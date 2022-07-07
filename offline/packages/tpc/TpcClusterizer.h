#ifndef TPC_TPCCLUSTERIZER_H
#define TPC_TPCCLUSTERIZER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsGeometry.h>

#include <map> 
#include <vector>
#include <string>

class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4CylinderCellGeom;
class PHG4CylinderCellGeomContainer;
class TrainingHitsContainer;

//typedef std::pair<int, int> iphiz;
//typedef std::pair<double, iphiz> ihit;
typedef std::pair<unsigned short, unsigned short> iphiz;
typedef std::pair<unsigned short, iphiz> ihit;

class TpcClusterizer : public SubsysReco
{
 public:
  TpcClusterizer(const std::string &name = "TpcClusterizer");
  ~TpcClusterizer() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_sector_fiducial_cut(const double cut){SectorFiducialCut = cut; }
  void set_use_nn(bool use_nn) { _use_nn = use_nn; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}
  void set_do_wedge_emulation(bool do_wedge){ do_wedge_emulation = do_wedge;}
  void set_max_cluster_half_size_phi(unsigned short size) { MaxClusterHalfSizePhi = size ;}
  void set_max_cluster_half_size_z(unsigned short size) { MaxClusterHalfSizeT = size ;}
  void set_cluster_version(int value) { cluster_version = value; }
  
 private:
  bool is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom) const;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  bool _use_nn = true;
  bool do_hit_assoc = true;
  bool do_wedge_emulation = false;
  double pedestal = 74.4;
  double SectorFiducialCut = 0.5;
  unsigned short MaxClusterHalfSizePhi = 3;
  unsigned short MaxClusterHalfSizeT = 5;
  int cluster_version = 3;
  double m_tdriftmax = 0;
  double AdcClockPeriod = 53.0;   // ns 

  // TPC shaping offset correction parameters
  // From Tony Frawley May 13, 2021
  double par0_neg = 0.0503;
  double par0_pos = -0.0503;
  
  TrainingHitsContainer *m_training;
};

#endif
