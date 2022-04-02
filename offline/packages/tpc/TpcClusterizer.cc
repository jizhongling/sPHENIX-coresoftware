#include "TpcClusterizer.h"

#include "TpcDefs.h"
#include "TrainingHitsContainer.h"
#include "TrainingHits.h"

#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                         // for PHIODataNode
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                             // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// One-stop header
#include <torch/script.h>

#include <TMatrixFfwd.h>    // for TMatrixF
#include <TMatrixT.h>       // for TMatrixT, ope...
#include <TMatrixTUtils.h>  // for TMatrixTRow

#include <TFile.h>  

#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <array>
#include <vector>
#include <algorithm>
// Terra incognita....
#include <pthread.h>

namespace 
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  
  typedef std::pair<unsigned short, unsigned short> iphiz;
  typedef std::pair<unsigned short, iphiz> ihit;
  typedef std::pair<TrkrDefs::cluskey, TrkrDefs::hitkey> assoc;
  
  // Neural network parameters and modules
  bool use_nn = true;
  const int nd = 5;
  const int nmax_truth = 5;
  torch::jit::script::Module module_ntruth;
  torch::jit::script::Module module_pos[nmax_truth];

  struct thread_data 
  {
    PHG4CylinderCellGeom *layergeom = nullptr;
    TrkrHitSet *hitset = nullptr;
    ActsSurfaceMaps *surfmaps = nullptr;
    ActsTrackingGeometry *tGeometry = nullptr;
    unsigned int layer = 0;
    int side = 0;
    unsigned int sector = 0;
    float pedestal = 0;
    bool do_assoc = true;
    unsigned short phibins = 0;
    unsigned short phioffset = 0;
    unsigned short zbins = 0;
    unsigned short zoffset = 0;
    unsigned short maxHalfSizeZ = 0;
    unsigned short maxHalfSizePhi = 0;
    double m_drift_velocity_scale = 1.0;
    double par0_neg = 0;
    double par0_pos = 0;
    std::map<TrkrDefs::cluskey, TrkrCluster *> *clusterlist = nullptr;
    //std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>  *clusterhitassoc = nullptr;
    std::vector<assoc> *association_vector = nullptr;
    std::vector<TrkrClusterv3*> *cluster_vector = nullptr;
    std::vector<TrainingHits*> *v_hits = nullptr;
  };
  
  pthread_mutex_t mythreadlock;
  
  void remove_hit(double adc, int phibin, int zbin, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
  {
    typedef std::multimap<unsigned short, ihit>::iterator hit_iterator;
    std::pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
    hit_iterator it = iterpair.first;
    for (; it != iterpair.second; ++it) {
      if (it->second.second.first == phibin && it->second.second.second == zbin) { 
	all_hit_map.erase(it);
	break;
      }
    }
    adcval[phibin][zbin] = 0;
  }
  
  void remove_hits(std::vector<ihit> &ihit_list, std::multimap<unsigned short, ihit> &all_hit_map,std::vector<std::vector<unsigned short>> &adcval)
  {
    for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
      unsigned short adc    = iter->first; 
      unsigned short phibin = iter->second.first;
      unsigned short zbin   = iter->second.second;
      remove_hit(adc,phibin,zbin,all_hit_map,adcval);
    }
  }
  
  void find_z_range(int phibin, int zbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, int& zdown, int& zup){
	
    const int FitRangeZ = (int) my_data.maxHalfSizeZ;
    const int NZBinsMax = (int) my_data.zbins;
	  zup = 0;
	  zdown = 0;
	  for(int iz=0; iz< FitRangeZ; iz++){
	    int cz = zbin + iz;
	
	    if(cz <= 0 || cz >= NZBinsMax){
	      // zup = iz;
	      break; // truncate edge
	    }
	    
	    if(adcval[phibin][cz] <= 0) {
	      break;
	    }
	    //check local minima and break at minimum.
	    if(cz<NZBinsMax-4){//make sure we stay clear from the edge
	      if(adcval[phibin][cz]+adcval[phibin][cz+1] < 
		 adcval[phibin][cz+2]+adcval[phibin][cz+3]){//rising again
		zup = iz+1;
		break;
	      }
	    }
	    zup = iz;
	  }
	  for(int iz=0; iz< FitRangeZ; iz++){
	    int cz = zbin - iz;
	    if(cz <= 0 || cz >= NZBinsMax){
	      //      zdown = iz;
	      break; // truncate edge
	    }
	    if(adcval[phibin][cz] <= 0) {
	      break;
	    }
	
	    if(cz>4){//make sure we stay clear from the edge
	      if(adcval[phibin][cz]+adcval[phibin][cz-1] < 
		 adcval[phibin][cz-2]+adcval[phibin][cz-3]){//rising again
		zdown = iz+1;
		break;
	      }
	    }
	    zdown = iz;
	  }
	}
	
  void find_phi_range(int phibin, int zbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, int& phidown, int& phiup)
  {
	
    int FitRangePHI = (int) my_data.maxHalfSizePhi;
    int NPhiBinsMax = (int) my_data.phibins;
	  phidown = 0;
	  phiup = 0;
	  for(int iphi=0; iphi< FitRangePHI; iphi++){
	    int cphi = phibin + iphi;
	    if(cphi < 0 || cphi >= NPhiBinsMax){
	      // phiup = iphi;
	      break; // truncate edge
	    }
	    
	    //break when below minimum
	    if(adcval[cphi][zbin] <= 0) {
	      // phiup = iphi;
	      break;
	    }
	    //check local minima and break at minimum.
	    if(cphi<NPhiBinsMax-4){//make sure we stay clear from the edge
	      if(adcval[cphi][zbin]+adcval[cphi+1][zbin] < 
		 adcval[cphi+2][zbin]+adcval[cphi+3][zbin]){//rising again
		phiup = iphi+1;
		break;
	      }
	    }
	    phiup = iphi;
	  }
	
	  for(int iphi=0; iphi< FitRangePHI; iphi++){
	    int cphi = phibin - iphi;
	    if(cphi < 0 || cphi >= NPhiBinsMax){
	      // phidown = iphi;
	      break; // truncate edge
	    }
	    
	    if(adcval[cphi][zbin] <= 0) {
	      //phidown = iphi;
	      break;
	    }
	
	    if(cphi>4){//make sure we stay clear from the edge
	      if(adcval[cphi][zbin]+adcval[cphi-1][zbin] < 
		 adcval[cphi-2][zbin]+adcval[cphi-3][zbin]){//rising again
		phidown = iphi+1;
		break;
	      }
	    }
	    phidown = iphi;
	  }
	}
	
  void get_cluster(int phibin, int zbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, std::vector<ihit> &ihit_list)
	{
          phibin = std::clamp(phibin, 0, (int)my_data.phibins - 1);
          zbin = std::clamp(zbin, 0, (int)my_data.zbins - 1);

	  // search along phi at the peak in z
	 
	  int zup =0;
	  int zdown =0;
	  find_z_range(phibin, zbin, my_data, adcval, zdown, zup);
	  //now we have the z extent of the cluster, go find the phi edges
	
	  for(int iz=zbin - zdown ; iz<= zbin + zup; iz++){
	    int phiup = 0;
	    int phidown = 0;
	    find_phi_range(phibin, iz, my_data, adcval, phidown, phiup);
	    for (int iphi = phibin - phidown; iphi <= (phibin + phiup); iphi++){
	      iphiz iCoord(std::make_pair(iphi,iz));
	      ihit  thisHit(adcval[iphi][iz],iCoord);
	      ihit_list.push_back(thisHit);
	    }
	  }
	}
  
  	Surface get_tpc_surface_from_coords(TrkrDefs::hitsetkey hitsetkey,
					    Acts::Vector3 world,
					    ActsSurfaceMaps *surfMaps,
					    ActsTrackingGeometry *tGeometry,
					    TrkrDefs::subsurfkey& subsurfkey)
	{

	  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
	  std::map<unsigned int, std::vector<Surface>>::iterator mapIter;
	  mapIter = surfMaps->tpcSurfaceMap.find(layer);
	  
	  if(mapIter == surfMaps->tpcSurfaceMap.end())
	    {
	      std::cout << PHWHERE 
			<< "Error: hitsetkey not found in clusterSurfaceMap, hitsetkey = "
			<< hitsetkey << std::endl;
	      return nullptr;
	    }
	
	  double world_phi = atan2(world[1], world[0]);
	  double world_z = world[2];
	  
	  std::vector<Surface> surf_vec = mapIter->second;
	  unsigned int surf_index = 999;

	  // Predict which surface index this phi and z will correspond to
	  // assumes that the vector elements are ordered positive z, -pi to pi, then negative z, -pi to pi
	  double fraction =  (world_phi + M_PI) / (2.0 * M_PI);
	  double rounded_nsurf = round( (double) (surf_vec.size()/2) * fraction  - 0.5);
	  unsigned int nsurf = (unsigned int) rounded_nsurf; 
	  if(world_z < 0)
	    nsurf += surf_vec.size()/2;

	  double surfStepPhi = tGeometry->tpcSurfStepPhi;
	  double surfStepZ = tGeometry->tpcSurfStepZ;

	  Surface this_surf = surf_vec[nsurf];
	  auto vec3d = this_surf->center(tGeometry->geoContext);
	  std::vector<double> surf_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
	  double surf_z = surf_center[2];
	  double surf_phi = atan2(surf_center[1], surf_center[0]);

	  if( ((world_phi >= (surf_phi - (surfStepPhi / 2.0))) && 
	       (world_phi <= (surf_phi + (surfStepPhi / 2.0)))) &&
	      ((world_z >= (surf_z - (surfStepZ / 2.0))) 
	       && (world_z <= (surf_z + (surfStepZ / 2.0)))) )	
	    {
	      surf_index = nsurf;
	      subsurfkey = nsurf;
	    }    
	  else
	    {
	      std::cout << PHWHERE 
			<< "Error: TPC surface index not defined, skipping cluster!" 
			<< std::endl;
	      std::cout << "     coordinates: " << world[0] << "  " << world[1] << "  " << world[2] 
			<< " radius " << sqrt(world[0]*world[0]+world[1]*world[1]) << std::endl;
	      std::cout << "     world_phi " << world_phi << " world_z " << world_z << std::endl;
	      std::cout << "     surf coords: " << surf_center[0] << "  " << surf_center[1] << "  " << surf_center[2] << std::endl;
	      std::cout << "     surf_phi " << surf_phi << " surf_z " << surf_z << std::endl; 
	      std::cout << "     surfStepPhi " << surfStepPhi << " surfStepZ " << surfStepZ << std::endl; 
	      std::cout << " number of surfaces " << surf_vec.size() << " nsurf: "  << nsurf << std::endl;
	      return nullptr;
	    }
	  	 
	  return surf_vec[surf_index];
	
	}
  
    void calc_cluster_parameter(const std::vector<ihit> &ihit_list,int iclus, const thread_data& my_data )
    {
    
      // get z range from layer geometry
      /* these are used for rescaling the drift velocity */
      const double z_min = -105.5;
      const double z_max = 105.5;
      
      // loop over the hits in this cluster
      double z_sum = 0.0;
      double phi_sum = 0.0;
      double adc_sum = 0.0;
      double z2_sum = 0.0;
      double phi2_sum = 0.0;
      
      double radius = my_data.layergeom->get_radius();  // returns center of layer
      
      int phibinhi = -1;
      int phibinlo = 666666;
      int zbinhi = -1;
      int zbinlo = 666666;
      int clus_size = ihit_list.size();
      
      if(clus_size == 1) return;
      
      std::vector<TrkrDefs::hitkey> hitkeyvec;
      for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
	double adc = iter->first; 
	
	if (adc <= 0) continue;
	
	int iphi = iter->second.first + my_data.phioffset;
	int iz   = iter->second.second + my_data.zoffset;
	if(iphi > phibinhi) phibinhi = iphi;
	if(iphi < phibinlo) phibinlo = iphi;
	if(iz > zbinhi) zbinhi = iz;
	if(iz < zbinlo) zbinlo = iz;
	
	// update phi sums
	double phi_center = my_data.layergeom->get_phicenter(iphi);
	phi_sum += phi_center * adc;
	phi2_sum += square(phi_center)*adc;
	
	// update z sums
	double z = my_data.layergeom->get_zcenter(iz);
	
	// apply drift velocity scale
	/* this formula ensures that z remains unchanged when located on one of the readout plane, at either z_max or z_min */
	z =
	  z*my_data.m_drift_velocity_scale +
	  (z>0 ? z_max:z_min)*(1.-my_data.m_drift_velocity_scale);
	
	z_sum += z*adc;
	z2_sum += square(z)*adc;
	
	adc_sum += adc;
	
	// capture the hitkeys for all adc values above a certain threshold
	TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
	// if(adc>5)
	hitkeyvec.push_back(hitkey);
      }
      if (adc_sum < 10){
	hitkeyvec.clear();
	return;  // skip obvious noise "clusters"
      }  
      // This is the global position
      double clusphi = phi_sum / adc_sum;
      double clusz = z_sum / adc_sum;
      float clusx = radius * cos(clusphi);
      float clusy = radius * sin(clusphi);
      
      const double phi_cov = phi2_sum/adc_sum - square(clusphi);
      const double z_cov = z2_sum/adc_sum - square(clusz);
       // Get the surface key to find the surface from the 

      TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey( my_data.layer, my_data.sector, my_data.side );      
      Acts::Vector3 global(clusx, clusy, clusz);
      TrkrDefs::subsurfkey subsurfkey = 0;
      Surface surface = get_tpc_surface_from_coords(tpcHitSetKey,
						    global,
						    my_data.surfmaps,
						    my_data.tGeometry,
						    subsurfkey);
      
      if(!surface)
	{
	  /// If the surface can't be found, we can't track with it. So 
	  /// just return and don't add the cluster to the container
	  hitkeyvec.clear();
	  return;
	}
      // create the cluster entry directly in the node tree
      
      const TrkrDefs::cluskey ckey = TpcDefs::genClusKey( my_data.hitset->getHitSetKey(), iclus );
      TrkrClusterv3 *clus = new TrkrClusterv3();
      //auto clus = std::make_unique<TrkrClusterv3>();
      clus->setClusKey(ckey);
      
      // Estimate the errors
      const double phi_err_square = (phibinhi == phibinlo) ?
	square(radius*my_data.layergeom->get_phistep())/12:
	square(radius)*phi_cov/(adc_sum*0.14);
      
      const double z_err_square = (zbinhi == zbinlo) ?
	square(my_data.layergeom->get_zstep())/12:
	z_cov/(adc_sum*0.14);
      
      // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
      // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
      //    - N is the number of electrons that drift to the readout plane
      // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
      // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
      // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
      // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14
      clusz -= (clusz<0) ? my_data.par0_neg:my_data.par0_pos;
      
      // Fill in the cluster details
      //================
      clus->setAdc(adc_sum);
      

      TMatrixF ERR(3, 3);
      ERR[0][0] = 0.0;
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = phi_err_square;  //cluster_v1 expects rad, arc, z as elementsof covariance
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = z_err_square;
      
     
      
      clus->setSubSurfKey(subsurfkey);
      
      Acts::Vector3 center = surface->center(my_data.tGeometry->geoContext)/Acts::UnitConstants::cm;
      
      // no conversion needed, only used in acts
      Acts::Vector3 normal = surface->normal(my_data.tGeometry->geoContext);
      double clusRadius = sqrt(clusx * clusx + clusy * clusy);
      double rClusPhi = clusRadius * clusphi;
      double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
      double surfPhiCenter = atan2(center[1], center[0]);
      double surfRphiCenter = surfPhiCenter * surfRadius;
      double surfZCenter = center[2];
      
      auto local = surface->globalToLocal(my_data.tGeometry->geoContext,
					  global * Acts::UnitConstants::cm,
					  normal);
      Acts::Vector2 localPos;
      
      // Prefer Acts transformation since we build the TPC surfaces manually
      if(local.ok())
	{
	  localPos = local.value() / Acts::UnitConstants::cm;
	}
      else
	{
	  /// otherwise take the manual calculation
	  localPos(0) = rClusPhi - surfRphiCenter;
	  localPos(1) = clusz - surfZCenter; 
	}
      clus->setLocalX(localPos(0));
      clus->setLocalY(localPos(1));
     
      // clus->setLocalX(clusx);
      // clus->setLocalY(clusy);
      clus->setActsLocalError(0,0, ERR[1][1]);
      clus->setActsLocalError(1,0, ERR[2][1]);
      clus->setActsLocalError(0,1, ERR[1][2]);
      clus->setActsLocalError(1,1, ERR[2][2]);

      // Add the hit associations to the TrkrClusterHitAssoc node
      // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)
      
      if(my_data.clusterlist)     
	{
	  my_data.cluster_vector->push_back(clus);
	  /*
	   * if cluster was properly inserted in the map, one should release the unique_ptr,
	   * to make sure that the cluster is not deleted when going out-of-scope
	   */
	  /*
	  if( inserted ) clus.release();
	  else {
	    // print error message. Duplicated cluster keys should not happen
	    std::cout 
	      << PHWHERE 
	      << "Error: duplicated cluster key: " << ckey << " - new cluster not inserted in map"
	      << std::endl;        
	  }
	  */
	  
	}

      //      if(my_data.do_assoc && my_data.clusterhitassoc){
      if(my_data.do_assoc){
	//std::cout << "filling assoc" << std::endl;
	for (unsigned int i = 0; i < hitkeyvec.size(); i++){
	  my_data.association_vector->push_back(std::make_pair(ckey, hitkeyvec[i]));
	}
      }
      hitkeyvec.clear();
    }
  
    void calc_cluster_parameter(float local_phi, float local_z,
        const std::vector<ihit> &ihit_list, int iclus, const thread_data& my_data)
    {
      local_phi = std::clamp(local_phi, 0.f, (float)my_data.phibins - 1.f);
      local_z = std::clamp(local_z, 0.f, (float)my_data.zbins - 1.f);
      int iphi = (int)local_phi + my_data.phioffset;
      int iz = (int)local_z + my_data.zoffset;

      // this is the global position
      double radius = my_data.layergeom->get_radius();  // return center of layer
      double dphi = (local_phi - std::floor(local_phi)) * my_data.layergeom->get_phistep();
      double dz = (local_z - std::floor(local_z)) * my_data.layergeom->get_zstep();
      double clusphi = my_data.layergeom->get_phicenter(iphi) + dphi;
      float clusx = radius * cos(clusphi);
      float clusy = radius * sin(clusphi);
      float clusz = my_data.layergeom->get_zcenter(iz) + dz;

      // get z range from layer geometry
      // these are used for rescaling the drift velocity
      const double z_min = -105.5;
      const double z_max = 105.5;

      // loop over the hits in this cluster
      double z_sum = 0.0;
      double phi_sum = 0.0;
      double adc_sum = 0.0;
      double z2_sum = 0.0;
      double phi2_sum = 0.0;

      int phibinhi = -1;
      int phibinlo = 666666;
      int zbinhi = -1;
      int zbinlo = 666666;
      int clus_size = ihit_list.size();

      if(clus_size <= 1) return;

      std::vector<TrkrDefs::hitkey> hitkeyvec;
      for(auto iter = ihit_list.begin(); iter != ihit_list.end(); iter++)
      {
	double adc = iter->first;
	if(adc <= 0) continue;

	int iphi = iter->second.first + my_data.phioffset;
	int iz   = iter->second.second + my_data.zoffset;
	if(iphi > phibinhi) phibinhi = iphi;
	if(iphi < phibinlo) phibinlo = iphi;
	if(iz > zbinhi) zbinhi = iz;
	if(iz < zbinlo) zbinlo = iz;

	// update phi sums
	double phi_center = my_data.layergeom->get_phicenter(iphi);
	phi_sum += phi_center * adc;
	phi2_sum += square(phi_center)*adc;

	// update z sums
	double z = my_data.layergeom->get_zcenter(iz);

	// apply drift velocity scale
	/* this formula ensures that z remains unchanged when located on one of the readout plane, at either z_max or z_min */
	z = z*my_data.m_drift_velocity_scale + (z>0 ? z_max:z_min)*(1.-my_data.m_drift_velocity_scale);

	z_sum += z*adc;
	z2_sum += square(z)*adc;

	adc_sum += adc;

	// capture the hitkeys for all adc values above a certain threshold
	TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
	hitkeyvec.push_back(hitkey);
      }

      if(adc_sum < 10)
      {
	hitkeyvec.clear();
	return;  // skip obvious noise "clusters"
      }

      // Get the surface key to find the surface
      TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(my_data.layer, my_data.sector, my_data.side);
      Acts::Vector3 global(clusx, clusy, clusz);
      TrkrDefs::subsurfkey subsurfkey = 0;
      Surface surface = get_tpc_surface_from_coords(tpcHitSetKey,
						    global,
						    my_data.surfmaps,
						    my_data.tGeometry,
						    subsurfkey);

      // If the surface can't be found, we can't track with it.
      // So just return and don't add the cluster to the container
      if(!surface) return;

      // Create the cluster entry directly in the node tree
      const TrkrDefs::cluskey ckey = TpcDefs::genClusKey(my_data.hitset->getHitSetKey(), iclus);
      TrkrClusterv3 *clus = new TrkrClusterv3();
      clus->setClusKey(ckey);

      // Estimate the errors
      // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
      // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
      //    - N is the number of electrons that drift to the readout plane
      // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
      // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
      // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
      // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14
      const double phi_cov = phi2_sum/adc_sum - square(phi_sum/adc_sum);
      const double z_cov = z2_sum/adc_sum - square(z_sum/adc_sum);

      const double phi_err_square = (phibinhi == phibinlo) ?
	square(radius*my_data.layergeom->get_phistep())/12 :
	square(radius)*phi_cov/(adc_sum*0.14);

      const double z_err_square = (zbinhi == zbinlo) ?
	square(my_data.layergeom->get_zstep())/12 :
	z_cov/(adc_sum*0.14);

      clusz -= (clusz < 0) ? my_data.par0_neg : my_data.par0_pos;

      // Fill in the cluster details
      clus->setAdc(adc_sum);
      clus->setSubSurfKey(subsurfkey);

      Acts::Vector3 center = surface->center(my_data.tGeometry->geoContext)/Acts::UnitConstants::cm;

      // No conversion needed, only used in acts
      Acts::Vector3 normal = surface->normal(my_data.tGeometry->geoContext);
      double clusRadius = sqrt(clusx * clusx + clusy * clusy);
      double rClusPhi = clusRadius * clusphi;
      double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
      double surfPhiCenter = atan2(center[1], center[0]);
      double surfRphiCenter = surfPhiCenter * surfRadius;
      double surfZCenter = center[2];

      auto local = surface->globalToLocal(my_data.tGeometry->geoContext,
					  global * Acts::UnitConstants::cm,
					  normal);
      Acts::Vector2 localPos;

      // Prefer Acts transformation since we build the TPC surfaces manually
      if(local.ok())
      {
        localPos = local.value() / Acts::UnitConstants::cm;
      }
      else
      {
        // Otherwise take the manual calculation
        localPos(0) = rClusPhi - surfRphiCenter;
        localPos(1) = clusz - surfZCenter;
      }
      clus->setLocalX(localPos(0));
      clus->setLocalY(localPos(1));

      clus->setActsLocalError(0,0, phi_err_square);
      clus->setActsLocalError(1,0, 0.);
      clus->setActsLocalError(0,1, 0.);
      clus->setActsLocalError(1,1, z_err_square);

      // Add the cluster to the cluster list
      if(my_data.clusterlist)
        my_data.cluster_vector->push_back(clus);

      // Add the hit associations to the TrkrClusterHitAssoc node
      // We need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)
      if(my_data.do_assoc)
	for (unsigned int i = 0; i < hitkeyvec.size(); i++)
	  my_data.association_vector->push_back(std::make_pair(ckey, hitkeyvec[i]));
      hitkeyvec.clear();
    }

  void *ProcessSector(void *threadarg) {

    auto my_data = static_cast<thread_data*>(threadarg);
    const auto& pedestal  = my_data->pedestal;
    const auto& phibins   = my_data->phibins;
    const auto& phioffset = my_data->phioffset;
    const auto& zbins     = my_data->zbins ;
    const auto& zoffset   = my_data->zoffset ;
    const auto& layer     = my_data->layer;
	
    TrkrHitSet *hitset = my_data->hitset;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    
    // for convenience, create a 2D vector to store adc values in and initialize to zero
    std::vector<std::vector<unsigned short>> adcval(phibins, std::vector<unsigned short>(zbins, 0));
    std::multimap<unsigned short, ihit> all_hit_map;
    std::vector<ihit> hit_vect;
    
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	 hitr != hitrangei.second;
	 ++hitr){
      if( TpcDefs::getPad(hitr->first) - phioffset < 0 ){
	//std::cout << "WARNING phibin out of range: " << TpcDefs::getPad(hitr->first) - phioffset << " | " << phibins << std::endl;
	continue;
      }
      if( TpcDefs::getTBin(hitr->first) - zoffset < 0 ){
	//std::cout << "WARNING zbin out of range: " << TpcDefs::getTBin(hitr->first) - zoffset  << " | " << zbins <<std::endl;
      }
      unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
      unsigned short zbin = TpcDefs::getTBin(hitr->first) - zoffset;
      if(phibin>=phibins){
	//std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
	continue;
      }
      if(zbin>=zbins){
	//std::cout << "WARNING z bin out of range: " << zbin << " | " << zbins << std::endl;
	continue;
      }
      float_t fadc = (hitr->second->getAdc()) - pedestal; // proper int rounding +0.5
      unsigned short adc = 0;
      if(fadc>0) adc =  (unsigned short) fadc;
      if(phibin >= phibins) continue;
      if(zbin   >= zbins) continue; // zbin is unsigned int, <0 cannot happen
      
      if(adc>0){
	if(adc>5){
	  iphiz iCoord(std::make_pair(phibin,zbin));
	  ihit  thisHit(adc,iCoord);
	  all_hit_map.insert(std::make_pair(adc, thisHit));
	}
	adcval[phibin][zbin] = (unsigned short) adc;
      }
    }
    
    std::vector<std::vector<unsigned short>> adcval_copy(adcval);
    std::multimap<unsigned short, ihit> all_hit_map_copy(all_hit_map);

    while(!use_nn && all_hit_map_copy.size() > 0)
    {
      auto training_hits = new TrainingHits;

      auto iter = all_hit_map_copy.rbegin();
      ihit hiHit = iter->second;
      int iphi = hiHit.second.first;
      int iz = hiHit.second.second;

      // This is the global position
      training_hits->phi = my_data->layergeom->get_phicenter(iphi+phioffset);
      training_hits->z = my_data->layergeom->get_zcenter(iz+zoffset);
      training_hits->z -= training_hits->z < 0 ? my_data->par0_neg : my_data->par0_pos;
      training_hits->layer = my_data->layer;

      for(int phibin = std::max(iphi-nd,0); phibin < std::min(iphi+nd+1,(int)phibins); phibin++)
        for(int zbin = std::max(iz-nd,0); zbin < std::min(iz+nd+1,(int)zbins); zbin++)
        {
          training_hits->v_adc[(phibin-iphi+nd)*(2*nd+1)+(zbin-iz+nd)] = adcval_copy[phibin][zbin];
          remove_hit(adcval_copy[phibin][zbin],phibin,zbin,all_hit_map_copy,adcval_copy);
        }

      my_data->v_hits->emplace_back(training_hits);
    }

    int nclus = 0;
    while(!use_nn && all_hit_map.size() > 0)
    {
      auto iter = all_hit_map.rbegin();
      ihit hiHit = iter->second;
      int iphi = hiHit.second.first;
      int iz = hiHit.second.second;

      // put all hits in the all_hit_map (sorted by adc)
      // start with highest adc hit
      // -> cluster around it and get vector of hits
      std::vector<ihit> ihit_list;
      get_cluster(iphi, iz, *my_data, adcval, ihit_list);

      // -> calculate cluster parameters
      // -> add hits to truth association
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      calc_cluster_parameter(ihit_list,nclus++, *my_data);
      remove_hits(ihit_list, all_hit_map, adcval);
      ihit_list.clear();
    }

    nclus = 0;
    while(use_nn && all_hit_map.size() > 0)
    {
      auto iter = all_hit_map.rbegin();
      ihit hiHit = iter->second;
      int iphi = hiHit.second.first;
      int iz = hiHit.second.second;
      int li = (layer - 7) / 16;
      float zr = my_data->layergeom->get_zcenter(iz+zoffset) / my_data->layergeom->get_radius();

      std::array<unsigned short, (2*nd+1)*(2*nd+1)> v_adc;
      v_adc.fill(0);
      for(int phibin = std::max(iphi-nd,0); phibin < std::min(iphi+nd+1,(int)phibins); phibin++)
        for(int zbin = std::max(iz-nd,0); zbin < std::min(iz+nd+1,(int)zbins); zbin++)
        {
          // put all hits in the 11x11 hits array (sorted by adc)
          // start with highest adc hit
          v_adc[(phibin-iphi+nd)*(2*nd+1)+(zbin-iz+nd)] = adcval[phibin][zbin];
          // remove hits from all_hit_map
          // repeat untill all_hit_map empty
          remove_hit(adcval[phibin][zbin],phibin,zbin,all_hit_map,adcval);
        }

      try
      {
        // Create a vector of inputs
        std::vector<torch::jit::IValue> inputs;
        inputs.emplace_back(torch::stack({
              torch::from_blob(std::vector<float>(v_adc.begin(), v_adc.end()).data(), {1, 2*nd+1, 2*nd+1}, torch::kFloat32),
              torch::full({1, 2*nd+1, 2*nd+1}, li, torch::kFloat32),
              torch::full({1, 2*nd+1, 2*nd+1}, zr, torch::kFloat32)
              }, 1));

        // Execute the model and turn its output into a tensor
        int ntruth = module_ntruth.forward(inputs).toTensor().argmax(1)[0].item<int>();
        if(ntruth < 1 || ntruth > nmax_truth) continue;
        at::Tensor ten_pos = module_pos[ntruth-1].forward(inputs).toTensor().clamp(-nd, nd);
        for(int i=0; i<ntruth; i++)
        {
          float clus_phi = ten_pos[0][0][i].item<float>() + static_cast<float>(iphi);
          float clus_z = ten_pos[0][1][i].item<float>() + static_cast<float>(iz);

          // -> cluster for the 11x11 hits array and get vector of hits
          // -> calculate cluster parameters
          // -> add hits to truth association
          std::vector<ihit> ihit_list;
          get_cluster(std::round(clus_phi), std::round(clus_z), *my_data, adcval_copy, ihit_list);
          calc_cluster_parameter(clus_phi, clus_z, ihit_list, nclus++, *my_data);
          ihit_list.clear();
        }
      }
      catch(const c10::Error &e)
      {
        std::cout << PHWHERE << "Error: Failed to execute NN modules" << std::endl;
      }
    }

    //std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey> empty_map;
    // empty_map.swap(*my_data.clusterhitassoc);
    // my_data->clusterhitassoc->clear();//HERE

    adcval.clear();
    all_hit_map.clear();
    adcval_copy.clear();
    all_hit_map_copy.clear();
    hit_vect.clear();
    pthread_exit(nullptr);
  }
}

TpcClusterizer::TpcClusterizer(const std::string &name)
  : SubsysReco(name)
{}

bool TpcClusterizer::is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom) const
{
  bool reject_it = false;

  // sector boundaries occur every 1/12 of the full phi bin range  
  int PhiBins = layergeom->get_phibins();
  int PhiBinsSector = PhiBins/12;

  double radius = layergeom->get_radius();
  double PhiBinSize = 2.0* radius * M_PI / (double) PhiBins;

  // sector starts where?
  int sector_lo = sector * PhiBinsSector;
  int sector_hi = sector_lo + PhiBinsSector - 1;

  int sector_fiducial_bins = (int) (SectorFiducialCut / PhiBinSize);

  if(phibin < sector_lo + sector_fiducial_bins || phibin > sector_hi - sector_fiducial_bins)
    {
      reject_it = true;
      /*
      int layer = layergeom->get_layer();
      std::cout << " local maximum is in sector fiducial boundary: layer " << layer << " radius " << radius << " sector " << sector 
      << " PhiBins " << PhiBins << " sector_fiducial_bins " << sector_fiducial_bins
      << " PhiBinSize " << PhiBinSize << " phibin " << phibin << " sector_lo " << sector_lo << " sector_hi " << sector_hi << std::endl;  
      */
    }

  return reject_it;
}


int TpcClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainerv3;
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  auto clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!clusterhitassoc)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    clusterhitassoc = new TrkrClusterHitAssocv3;
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  auto training_container = findNode::getClass<TrainingHitsContainer>(dstNode, "TRAINING_HITSET");
  if (!training_container)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    training_container = new TrainingHitsContainer;
    PHIODataNode<PHObject> *TrainingHitsContainerNode =
        new PHIODataNode<PHObject>(training_container, "TRAINING_HITSET", "PHObject");
    DetNode->addNode(TrainingHitsContainerNode);
  }

  use_nn = _use_nn;
  if(use_nn)
  {
    const char *fileroot = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/zji/github/sphenix/coresoftware/offline/packages/tpc/modules/net_model";
    char filename[200];
    try
    {
      // Deserialize the ScriptModule from a file using torch::jit::load()
      sprintf(filename, "%s-type0-nout%d.pt", fileroot, nmax_truth);
      module_ntruth = torch::jit::load(filename);
      for(int ntruth=1; ntruth<=nmax_truth; ntruth++)
      {
        sprintf(filename, "%s-type1-nout%d.pt", fileroot, ntruth);
        module_pos[ntruth-1] = torch::jit::load(filename);
      }
      std::cout << PHWHERE << "Use NN clustering" << std::endl;
    }
    catch(const c10::Error &e)
    {
      std::cout << PHWHERE << "Error: Cannot load module " << filename << std::endl;
      exit(1);
    }
  }
  else
  {
    std::cout << PHWHERE << "Use traditional clustering" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::process_event(PHCompositeNode *topNode)
{
  //  int print_layer = 18;

  if (Verbosity() > 1000)
    std::cout << "TpcClusterizer::Process_Event" << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainerv3>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssocv3>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for training hits
  m_training = findNode::getClass<TrainingHitsContainer>(topNode, "TRAINING_HITSET");
  if (!m_training)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRAINING_HITSET." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4CylinderCellGeomContainer *geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,
							 "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE
		<< "ActsTrackingGeometry not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode,
						   "ActsSurfaceMaps");
  if(!m_surfMaps)
    {
      std::cout << PHWHERE 
		<< "ActsSurfaceMaps not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters

  // loop over the TPC HitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  const int num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);

  // create structure to store given thread and associated data
  struct thread_pair_t
  {
    pthread_t thread;
    thread_data data;
  };
  
  // create vector of thread pairs and reserve the right size upfront to avoid reallocation
  std::vector<thread_pair_t> threads;
  threads.reserve( num_hitsets );

  std::vector<assoc> association_vectors;

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  if (pthread_mutex_init(&mythreadlock, nullptr) != 0)
    {
      printf("\n mutex init failed\n");
      return 1;
    }
  int count = 0;
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    //if(count>0)continue;
    const auto hitsetid = hitsetitr->first;
    // std::cout << " starting thread # " << count << std::endl;
    TrkrHitSet *hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
    PHG4CylinderCellGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    
    // instanciate new thread pair, at the end of thread vector
    thread_pair_t& thread_pair = threads.emplace_back();
    
    thread_pair.data.layergeom = layergeom;
    thread_pair.data.hitset = hitset;
    thread_pair.data.layer = layer;
    thread_pair.data.pedestal = pedestal;
    thread_pair.data.sector = sector;
    thread_pair.data.side = side;
    thread_pair.data.do_assoc = do_hit_assoc;
    thread_pair.data.clusterlist = m_clusterlist->getClusterMap(hitsetid);
    // thread_pair.data.clusterhitassoc = m_clusterhitassoc->getClusterMap(hitsetid);
    // thread_pair.data.clusterhitassoc = m_clusterhitassoc;
    thread_pair.data.association_vector  = new std::vector<assoc>;
    thread_pair.data.cluster_vector  = new std::vector<TrkrClusterv3*>;
    thread_pair.data.v_hits = new std::vector<TrainingHits*>;
    thread_pair.data.tGeometry = m_tGeometry;
    thread_pair.data.surfmaps = m_surfMaps;
    thread_pair.data.maxHalfSizeZ =  MaxClusterHalfSizeZ;
    thread_pair.data.maxHalfSizePhi = MaxClusterHalfSizePhi;
    thread_pair.data.m_drift_velocity_scale = m_drift_velocity_scale;
    thread_pair.data.par0_neg = par0_neg;
    thread_pair.data.par0_pos = par0_pos;
    
    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins/12;
    unsigned short NZBins = (unsigned short)layergeom->get_zbins();
    unsigned short NZBinsSide = NZBins/2;
    unsigned short NZBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;

    if (side == 0){
      NZBinsMin = 0;
    }
    else{
      NZBinsMin = NZBins / 2;
    }

    unsigned short ZOffset = NZBinsMin;

    thread_pair.data.phibins   = NPhiBinsSector;
    thread_pair.data.phioffset = PhiOffset;
    thread_pair.data.zbins     = NZBinsSide;
    thread_pair.data.zoffset   = ZOffset ;

    int rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *)&thread_pair.data);
    if (rc) {
      std::cout << "Error:unable to create thread," << rc << std::endl;
    }
    count++;
  }
  
  pthread_attr_destroy(&attr);
  count =0;
  // wait for completion of all threads
  for( const auto& thread_pair:threads )
  { 
    int rc2 = pthread_join(thread_pair.thread, nullptr);
    if (rc2) 
    { std::cout << "Error:unable to join," << rc2 << std::endl; }
    for(auto iter = thread_pair.data.association_vector->begin(); iter != thread_pair.data.association_vector->end();++iter){
      TrkrDefs::cluskey ckey = iter->first;
      TrkrDefs::hitkey hkey = iter->second;
      m_clusterhitassoc->addAssoc(ckey,hkey);
    }
    thread_pair.data.association_vector->clear();
    delete thread_pair.data.association_vector;

    for(auto citer = thread_pair.data.cluster_vector->begin(); citer != thread_pair.data.cluster_vector->end();++citer){
      m_clusterlist->addCluster(*citer);
      // m_clusterlist->addClusToVec(*citer);
      
      // delete *citer;
      //      const auto [iter, inserted] = my_data.clusterlist->insert(std::make_pair(ckey, clus.get()));
	  
	  /*
	   * if cluster was properly inserted in the map, one should release the unique_ptr,
	   * to make sure that the cluster is not deleted when going out-of-scope
	   */
    }
    thread_pair.data.cluster_vector->clear();
    delete thread_pair.data.cluster_vector;

    for(auto hiter = thread_pair.data.v_hits->begin(); hiter != thread_pair.data.v_hits->end(); hiter++)
      m_training->v_hits.emplace_back(**hiter);
    thread_pair.data.v_hits->clear();
    delete thread_pair.data.v_hits;
  }

  if (Verbosity() > 0)
    std::cout << "TPC Clusterizer found " << m_clusterlist->size() << " Clusters "  << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
