#include "TruthRecoTrackMatching.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>  // for PHDataNode
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <trackbase/EmbRecoMatch.h>
#include <trackbase/EmbRecoMatchContainer.h>
#include <trackbase/EmbRecoMatchContainerv1.h>
#include <trackbase/EmbRecoMatchv1.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrTruthTrack.h>
#include <trackbase/TrkrTruthTrackContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TSystem.h>

#include <algorithm>

// To change:
// Make the definition of matching clusters to be if the truth cluster center is withing 1/2 of width of the reco track center

typedef TruthRecoTrackMatching::RECO_pair_iter RECO_pair_iter;
typedef TruthRecoTrackMatching::PossibleMatch  PossibleMatch;

using std::make_tuple;
using std::make_pair;
using std::vector;
using std::get;
using std::array;
using std::endl;
using std::cout;

//! square
template<class T> inline constexpr T square( const T& x ) { return x*x; }

//! radius
template<class T> inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

//! pt
template<class T> T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

//! p
template<class T> T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

//! eta
template<class T> T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

TruthRecoTrackMatching::TruthRecoTrackMatching
      ( const unsigned short _nmincluster_match
      , const float _nmincluster_ratio
      , const double _cutoff_dphi
      , const double _same_dphi
      , const double _cutoff_deta
      , const double _same_deta
      , const double _cluster_nzwidths
      , const double _cluster_nphiwidths
      , const unsigned short    _max_nreco_per_truth
      , const unsigned short    _max_ntruth_per_reco
      ) : m_nmincluster_match { _nmincluster_match  } // minimum number of clusters to match, default=4
        , m_nmincluster_ratio { _nmincluster_ratio } // minimum ratio to match a track, default=0.
          // -- Track Kinematic Cuts to match --
        , m_cutoff_dphi         { _cutoff_dphi         } // maximum |dphi|=|phi_reco-phi_truth| to search for match
        , m_same_dphi           { _same_dphi           } // all tracks in this |dphi| must be tested for matches
        , m_cutoff_deta         { _cutoff_deta         } // same as m_cutoff_dphi for deta
        , m_same_deta           { _same_deta           } // same as m_same_dphi for deta
          // cluster matching widths (how close the truth center must be reco center)          
        , m_cluster_nzwidths    { _cluster_nzwidths    }
        , m_cluster_nphiwidths  { _cluster_nphiwidths  }
          // number of truth tracks allowed matched per reco track, and v. versa
        , m_max_nreco_per_truth { _max_nreco_per_truth }
        , m_max_ntruth_per_reco { _max_ntruth_per_reco }
        , m_nmatched_index_true {}
        , m_nmatched_id_reco {nullptr}
        , m_nmatched_id_true {nullptr}
{ 
  if (Verbosity() > 50) cout << " Starting TruthRecoTrackMatching.cc " << endl;
}

int TruthRecoTrackMatching::InitRun(PHCompositeNode *topNode) //`
{ 
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_nmatched_id_reco = &(m_EmbRecoMatchContainer->map_nTruthPerReco());
  m_nmatched_id_true = &(m_EmbRecoMatchContainer->map_nRecoPerTruth());
  return Fun4AllReturnCodes::EVENT_OK;
  if (Verbosity()>50) topNode->print();
  return 0;
}


int TruthRecoTrackMatching::process_event(PHCompositeNode* topnode)  //`
{
  if (topnode == nullptr) {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (Verbosity()>1000) topnode->print(); // perhaps not needed

  // cleanup output container
  if (m_TrackEvalContainer) {
    m_TrackEvalContainer->Reset();
    m_TrackEvalContainer->clearEvents();
    m_TrackEvalContainer->clearClusters();
    m_TrackEvalContainer->clearTracks();
  }

  // -------------------------------------------------------------------------------
  // Build recoData
  // ------------------------------------------------------------------------------

  // recoData is a vector of tuples that acts like a table with four columns, 
  // with one entry each of n tracks:
  //    (0)::float  (1)::float  (2)::float  (3)::short
  //     phi-0        eta-0       pT-0       index-0
  //     phi-1        eta-1       pT-1       index-1
  //     phi-2        eta-2       pT-2       index-2
  //     ...          ...         ...        ...   
  //     phi-n-2      eta-n-2     pT-n-2     index-n-2
  //     phi-n-1      eta-n-1     pT-n-1     index-n-1
  //     phi-n        eta-n       pT-n       index-n
  if (Verbosity()>60) std::cout << "reco tracks size: " << m_SvtxTrackMap->size() << std::endl;

  recoData.clear();


  for (auto reco = m_SvtxTrackMap->begin(); reco != m_SvtxTrackMap->end(); ++reco) {
    auto index_reco = reco->first;
    auto track      = reco->second;
    recoData.push_back( {track->get_phi(), track->get_eta(), track->get_pt(), index_reco } );
  }
  // sort the recoData table by phi
  std::sort(recoData.begin(), recoData.end(), CompRECOtoPhi());

  // phi will be sorted by proximity in the table, so re-add the first entries to the
  // end of the table so that they can be compared against the phi across the 
  // circular boundary condition. This makes the table potentially look like:
  //    (0)::float  (1)::float  (2)::float  (3)::short
  //     phi-0        eta-0       pT-0       index-0
  //     phi-1        eta-1       pT-1       index-1
  //     phi-2        eta-2       pT-2       index-2
  //     ...          ...         ...        ...   
  //     phi-n-2      eta-n-2     pT-n-2     index-n-2
  //     phi-n-1      eta-n-1     pT-n-1     index-n-1
  //     phi-n        eta-n       pT-n       index-n
  //     phi-0        eta-0       pT-0       index-0 // <- values wrapped from the lowest phi values
  //     phi-1        eta-1       pT-1       index-1 // <-
  RECOvec wraps{};
  auto wrap_from_start = std::upper_bound(recoData.begin(), recoData.end(), (-M_PI+m_cutoff_dphi),  CompRECOtoPhi());
  for (auto iter = recoData.begin(); iter != wrap_from_start; ++iter) {
    auto entry = *iter; // make a new copy to wrap to the other side of recoData
    get<RECOphi>(entry) = get<RECOphi>(entry) + 2*M_PI; // put the new copy on the other end
    wraps.push_back(entry);
  }
  for (auto E : wraps) recoData.push_back(E);

  if (Verbosity() > 70) {
    cout << " ****************************************** " << endl;
    cout << " This is the RECO map of tracks to match to " << endl;
    cout << " ****************************************** " << endl;
    for (auto& E : recoData) { 
      cout << Form(" id:%2i  (phi:eta:pt) (%5.2f:%5.2f:%5.2f)", get<RECOid>(E),
          get<RECOphi>(E), get<RECOeta>(E), get<RECOpt>(E)) << endl; }
    cout << " ****end*listing*********************** " << endl;
  }

  /******************************************************************************
   * Loop through the truth tracks one at a time. Based on track phi, eta, and pT
   * build two indices of truth-to-reco pairs:
   *    innerbox_pairs: pairs of truth-to-reco tracks within same_d{phi,eta}
   *      i.e. |phi_true-phi_reco|<same_dphi && |eta_true-eta_reco|<same_deta).
   *      These are the "innerboxes" (retangles in phi and eta space).
   *      All of these pairs will have to be checked for matching tracks, and the
   *      best fits will be removed first.
   *    outerbox_pairs: these are wider boxes (sized by cutoff_d{phi,eta})
   *      These are possible matches that are only checked for tracks remaining
   *      after the innerbox_pairs are checked and matches made.
   ******************************************************************************/

  std::vector<std::pair<unsigned short, unsigned short>> outerbox_pairs {};
  std::vector<std::pair<unsigned short, unsigned short>> innerbox_pairs {};

  if (Verbosity()>70) std::cout << "Number of truth tracks: " << m_TrkrTruthTrackContainer->getTruthTracks().size() << std::endl;

  unsigned short index_true {0};
  for (auto track : m_TrkrTruthTrackContainer->getTruthTracks()) {
    auto match_indices = find_box_matches(track->getPhi(), track->getPseudoRapidity(), track->getPt()); 

    // keep track of all truth tracks (by track-id) which has been matched
    for (auto& id_reco : match_indices.first)  innerbox_pairs.push_back( {index_true, id_reco} );
    for (auto& id_reco : match_indices.second) outerbox_pairs.push_back( {index_true, id_reco} );
    ++index_true;
  }

  match_tracks_in_box (innerbox_pairs);
  match_tracks_in_box (outerbox_pairs);


  // loop through all the truth track id's in order to fill the EmbRecoMatchContainer with the list of un-matched tracks
  for (auto track : m_TrkrTruthTrackContainer->getTruthTracks()) {
    m_EmbRecoMatchContainer->checkfill_idsTruthUnmatched(track->getTrackid());
  }

  if (Verbosity()>90) {
    cout << " --0-- Printing all matches stored (start)" << endl;
  // re-print all the tracks with the matches with the fit values
    for (auto match : m_EmbRecoMatchContainer->getMatches()) {
      cout << Form(" Match id(%2i->%2i) nClusMatch-nClusTrue-nClusReco (%2i:%2i:%2i)",
          match->idTruthTrack(), match->idRecoTrack(),
          match->nClustersMatched(), match->nClustersTruth(), match->nClustersReco()) << endl;
    }
    cout << " --1-- Printing all matches stored (done)" << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthRecoTrackMatching::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

  //--------------------------------------------------
  // Internal functions
  //--------------------------------------------------

int TruthRecoTrackMatching::createNodes(PHCompositeNode* topNode)
{
  // Initiailize the modules
  m_PHG4TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_PHG4TruthInfoContainer)
  {
    std::cout << "Could not locate G4TruthInfo node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_SvtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_SvtxTrackMap)
  {
    std::cout << "Could not locate SvtxTrackMap node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  m_TruthClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_TruthClusterContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHCLUSTERCONTAINER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_RecoClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_RecoClusterContainer)
  {
    std::cout << "Could not locate TRKR_CLUSTER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_TrkrTruthTrackContainer = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!m_TrkrTruthTrackContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHTRACKCONTAINER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_PHG4TpcCylinderGeomContainer = 
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_PHG4TpcCylinderGeomContainer)
  {
    std::cout << "Could not locate CYLINDERCELLGEOM_SVTX node when running \"TruthRecoTrackMatching\" module." << std::endl;
    /* std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl; */
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // note that layers 0-6, and > 55, don't work
  for (int layer=7; layer<55; ++layer) {
    PHG4TpcCylinderGeom *layergeom = m_PHG4TpcCylinderGeomContainer->GetLayerCellGeom(layer);
    if (layer==7) m_zstep = layergeom->get_zstep();
    m_phistep[layer] = layergeom->get_phistep();
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode,"TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer) {
    PHNodeIterator dstiter(dstNode);

    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    m_EmbRecoMatchContainer = new EmbRecoMatchContainerv1;
    auto newNode = new PHIODataNode<PHObject>(m_EmbRecoMatchContainer, "TRKR_EMBRECOMATCHCONTAINER", "PHObject");
    DetNode->addNode(newNode);

    m_TrackEvalContainer = new TrackEvaluationContainerv1;
    auto newNode2 = new PHIODataNode<PHObject>(m_TrackEvalContainer, "TrackMatchEvalContainer", "PHObject");
    DetNode->addNode(newNode2);
  }
  
  m_ActsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  return Fun4AllReturnCodes::EVENT_OK;
}

std::pair<std::vector<unsigned short>, std::vector<unsigned short>> 
TruthRecoTrackMatching::find_box_matches(float truth_phi, float truth_eta, float truth_pt) {
  // sort through the recoData to find:
  //     inner_box : possible matches withing m_same_d{phi,eta,pt}
  //     outer_box : (assumed to be strictly larger than, and therefore contain, inner_box, in m_cutoff_d{phi, eta,pt}
  // recoData is already sorted by eta, but as the boxes "zoom in" on acceptable matches, it gets
  // sorted by eta and (if applicable) pT as well. `mix_pair` keeps track of the range of recoData that has
  // been sorted so that it can be resorted to phi order.
  RECO_pair_iter mix_pair { recoData.begin(), recoData.end() };

  mix_pair.first = std::lower_bound(mix_pair.first, mix_pair.second, 
      truth_phi-m_cutoff_dphi, CompRECOtoPhi());

  mix_pair.second = std::upper_bound(mix_pair.first, mix_pair.second,
      truth_phi+m_cutoff_dphi, CompRECOtoPhi());

  if (mix_pair.first == mix_pair.second) {
    // there are no possible phi_matches; return nothing
    return { {}, {} };
  }
  
  // There are tracks within the outer_box phi range;
  // now re-shuffle the phi outerbox range for eta and see if there
  // are tracks in outer_box-eta range
  std::sort(mix_pair.first, mix_pair.second, CompRECOtoEta());

  RECO_pair_iter outer_box = mix_pair;
  outer_box.first = std::lower_bound(mix_pair.first, mix_pair.second,
      truth_eta-m_cutoff_deta, CompRECOtoEta());

  outer_box.second = std::lower_bound(outer_box.first, outer_box.second,
      truth_eta+m_cutoff_deta, CompRECOtoEta());

  if (outer_box.first == outer_box.second) {
    // there are no eta_matches in outerbox; resort mix_pair and return nothing
    std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());
    return { {}, {} };
  }

  // if pt is specified, restrict pt_box to the given range
  RECO_pair_iter inner_box = outer_box;
  const float _delta_outer_pt = delta_outer_pt(truth_pt);
  const float _delta_inner_pt = delta_inner_pt(truth_pt);
  if (_delta_outer_pt > 0) {
    std::sort(outer_box.first, outer_box.second, CompRECOtoPt());
    outer_box.first  = std::lower_bound(outer_box.first, outer_box.second, truth_pt-_delta_outer_pt, CompRECOtoPt());
    outer_box.second = std::upper_bound(outer_box.first, outer_box.second, truth_pt+_delta_outer_pt, CompRECOtoPt());

    // if not outer_box, resort and return nothing
    if (outer_box.first == outer_box.second) {
      // there are no eta_matches in outerbox; resort mix_pair and return nothing
      std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());
      return { {}, {} };
    }

    // now calculate the inner box -- this time sorting in reverse -- pT, then eta, then phi
    inner_box = outer_box;
    if (_delta_inner_pt > 0) {
      // start the inner pair -- the outerbox is already sorted by pT
      inner_box.first  = std::lower_bound(inner_box.first, inner_box.second, truth_pt-_delta_inner_pt, CompRECOtoPt());
      inner_box.second = std::upper_bound(inner_box.first, inner_box.second, truth_pt+_delta_inner_pt, CompRECOtoPt());
    }

    // go back to sorted by eta
    std::sort(inner_box.first, inner_box.second, CompRECOtoEta());
  }

  // At this point we know that we have outer_box matches
  // Finish checking if there are any possible inner_box matches

  // check for inner_box eta matches
  if (inner_box.first != inner_box.second) {
    inner_box.first = std::lower_bound(inner_box.first, inner_box.second,
        truth_eta-m_same_deta, CompRECOtoEta());

    inner_box.second = std::lower_bound(inner_box.first, inner_box.second,
        truth_eta+m_same_deta, CompRECOtoEta());
  }

  // check for inner_box phi matches
  if (inner_box.first != inner_box.second) {
    std::sort(inner_box.first, inner_box.second, CompRECOtoPhi());

    inner_box.first = std::lower_bound(inner_box.first, inner_box.second,
        truth_phi-m_same_dphi, CompRECOtoPhi());

    inner_box.second = std::lower_bound(inner_box.first, inner_box.second,
        truth_phi+m_same_dphi, CompRECOtoPhi());
  }

  // At this point there are definitely outer_box matches, and
  // possible inner_box matches. 
  // Return all these possible matches, evaluating all inner_box_matches
  std::vector<unsigned short> inner_box_matches, outer_box_matches;

  // fill inner_box_matches
  for (auto iter = inner_box.first; iter != inner_box.second; ++iter) {
    inner_box_matches.push_back( get<RECOid>(*iter) );
  }
  // fill inner_box_matches
  for (auto iter = outer_box.first; iter != inner_box.first; ++iter) {
    outer_box_matches.push_back(get<RECOid>(*iter));
  }
  for (auto iter = inner_box.second; iter != outer_box.second; ++iter) {
    outer_box_matches.push_back(get<RECOid>(*iter));
  }
  
  // resort recoData back to phi order
  std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());

  // return the box matches
  return std::make_pair(inner_box_matches, outer_box_matches);
}

void TruthRecoTrackMatching::match_tracks_in_box(
      std::vector<std::pair<unsigned short,unsigned short>>& box_pairs // possible matches by proximity in eta, phi, pT
) {
  if (box_pairs.size() == 0) return;

  std::sort(box_pairs.begin(), box_pairs.end()); // sorted by first index_true, then id_reco

  vector<PossibleMatch> poss_matches;
  vector<std::map<TrkrDefs::cluskey,TrkrDefs::cluskey>> v_clusmap; // key_reco, key_true
  unsigned short index_clusmap = 0;

  auto ipair = box_pairs.begin(); // keep track of current examined pair
  while (ipair != box_pairs.end()) {
    auto index_true = ipair->first;

    // See if this track already has the maximum number of matches. If it does, then skip all pairs using this truth-track-index;
    if (at_nmax_index_true(index_true)) {
      while (ipair != box_pairs.end()) {
        ++ipair;
        if (ipair->first != index_true) break;
      }
      continue;
    }

    // make all possible reco matches matched for this track
    std::map<TrkrDefs::hitsetkey,TrkrDefs::cluskey> truth_keys;
    auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[index_true];
    for (auto& key : truth_track->getClusters()) {
      truth_keys[TrkrDefs::getHitSetKeyFromClusKey(key)] = key;
    }
    unsigned short nclus_true = truth_keys.size();

    while (ipair != box_pairs.end()) { // Loop over all possible matches only for this index_true (subsequent true tracks will be caught on following loops)
      if (ipair->first != index_true) break; // never true on first iteration
      if (at_nmax_id_reco(ipair->second)) {
        ++ipair;
        continue;
      }// make sure the reco-track isn't alread matched
      // ok, make a possible match: compare the clusters in the truth track and the reco track
      unsigned short nclus_match   = 0; // fill in the comparison loop
      unsigned short nclus_nomatch = 0; // number of reco and truth cluster that share a hitsetkey, but still fair matching criteria
      unsigned short nclus_reco    = 0; // count in the comparison loop
      std::map<TrkrDefs::cluskey,TrkrDefs::cluskey> clus_map; // key_reco, key_true
      SvtxTrack* reco_track = m_SvtxTrackMap->get(ipair->second);
      auto tpcseed = reco_track->get_tpc_seed();
      if (tpcseed) {
        for (auto r_key = tpcseed->begin_cluster_keys(); r_key!=tpcseed->end_cluster_keys(); ++r_key) {
          ++nclus_reco;
          auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(*r_key);
          if (truth_keys.count(hitsetkey) != 0) {
            // reco and truth cluster are in same hitsetkey-indexed subsurface. See if they match (++nclus_match) or not (++nclus_nomatch)
            if (compare_cluster_pair(truth_keys[hitsetkey], *r_key, hitsetkey).first) {
              ++nclus_match;
              if (clus_map.count(*r_key) == 0)
                clus_map.insert({*r_key, truth_keys[hitsetkey]});
            } else {
              ++nclus_nomatch;
            }
          }
        }
      } // end reco seed loop
      // if the match passes minimum cuts, then it is a possible match
      if (Verbosity()>100) {
        auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[index_true];
        cout << Form(" Poss. pair: (phi,eta,pT) Match true(%7.4f,%7.4f,%7.4f) reco(%7.4f,%7.4f,%7.4f) delta(%7.4f,%7.4f,%7.4f) "
            "(nClMatch:nClTrue:nClReco:nClnomatch)(%2i-%2i-%2i-%2i)", 
            truth_track->getPhi(), truth_track->getPseudoRapidity(), truth_track->getPt(),
            reco_track->get_phi(), reco_track->get_eta(),            reco_track->get_pt(),
            abs_dphi(truth_track->getPhi(),reco_track->get_phi()),
            abs(truth_track->getPseudoRapidity()-reco_track->get_eta()),
            abs(truth_track->getPt() -reco_track->get_pt()),
            (int)nclus_match, (int)nclus_true, (int)nclus_reco, (int)nclus_nomatch);
      }
      if ( nclus_match >= m_nmincluster_match 
          && ( static_cast<float>(nclus_match)/nclus_true >= m_nmincluster_ratio)
         ) {
        if (Verbosity()>100) cout << " Y" << endl;
        poss_matches.push_back( {nclus_match, nclus_true, nclus_reco, ipair->first, ipair->second, index_clusmap} );
        v_clusmap.push_back(clus_map);
        index_clusmap++;
      } else {
        if (Verbosity()>100) cout << " N" << endl;
      }
      ++ipair;
    }
  }

  // ok add all possible matches started for the largest PM_nmatch (the top)
  //  for groups of ties of PM_nmatch, sort by PM_ntrue (from smallest)
  //    for groups of ties (PM_nmatch, PM_ntrue), go by smallest PM_nreco
  //       for groups of ties (PM_nmatch, PM_ntrue, PM_nreco), do a detailed sum diff in the deltas
  std::sort(poss_matches.begin(), poss_matches.end(), SortPossibleMatch());

  if (Verbosity()>200) {
    cout << " All Y possible matches (" << poss_matches.size() << ") track pairs  (nClMatched-nClTrue-nClReco : idTrue-idReco)" << endl;
    int i{0};
    for (auto match : poss_matches) {
      /* auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[get<PM_idtrue>(match)]; */
      /* auto index_trut = truth_track->getTrackid(); */
      cout << Form(" pair(%2i):  %2i-%2i-%2i-<%2i>-%2i ", i++
      , get<PM_nmatch> (match) 
      , get<PM_ntrue>  (match) 
      , get<PM_nreco>  (match) 
      , get<PM_idtrue> (match)
      , get<PM_idreco> (match) ) << endl;
    }
  }

  /* std::set<int> matched_idreco, matched_idtrue; */
  auto iter = poss_matches.begin();
  while (iter != poss_matches.end()) {
    if ( skip_match(*iter)) { 
      ++iter;
      continue;
    }
    vector<std::pair<float,int>> sigma_metric = {{0.,0}};
    int n_sigma = 0;
    auto iter_same = iter+1; // iterator to look forward from first point and get detailed comparisons for all possibly matched tracks
    while (
          iter_same != poss_matches.end()
       && (*iter_same) [PM_nmatch] == (*iter)[PM_nmatch]
       && (*iter_same) [PM_ntrue]  == (*iter)[PM_ntrue]
       && (*iter_same) [PM_nreco]  == (*iter)[PM_nreco]
    ) {
      ++n_sigma;
      if (n_sigma==1) { sigma_metric[0].first = sigma_CompMatchClusters(*iter); }
      sigma_metric.push_back({sigma_CompMatchClusters(*iter_same), n_sigma});
      ++iter_same;
    }
    std::sort(sigma_metric.begin(), sigma_metric.end());

    bool first = true;
    for (auto& sigma : sigma_metric) {
      if (first) first = false;
      else if (skip_match(*(iter+sigma.second))) continue;
      auto match = *(iter+sigma.second);
      auto index_true = match[PM_idtrue];
      auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[index_true];
      auto id_reco = match[PM_idreco];
      auto id_true = truth_track->getTrackid();
      auto save_match = new EmbRecoMatchv1( id_true, id_reco,
          match[PM_ntrue], match[PM_nreco], match[PM_nmatch]);
      m_EmbRecoMatchContainer->addMatch(save_match);
      add_match_eval(id_reco, id_true, v_clusmap[match[PM_clusmap]]);

      if (m_nmatched_index_true.find(index_true) == m_nmatched_index_true.end()) {
        m_nmatched_index_true[index_true] = 1;
      } else {
        m_nmatched_index_true[index_true] += 1;
      }

    }
    iter += sigma_metric.size();
  }
}

// ----------------------------------------
// convenience function
// ----------------------------------------
inline float TruthRecoTrackMatching::abs_dphi (float phi0, float phi1) {
  float dphi = fabs(phi0-phi1);
  while (dphi > M_PI) dphi -= fabs(2*M_PI);
  return dphi;
}

/* 
  input: TrkrCluster* truth, TrkrCluster *reco, TrkrDefs::hitsetkey, bool calc_sigma=false
  pair<bool, float> 
  out:   Does it match? 

 */
std::pair<bool, float> TruthRecoTrackMatching::compare_cluster_pair (
    TrkrDefs::cluskey key_T
  , TrkrDefs::cluskey key_R
  , TrkrDefs::hitsetkey hitsetkey
  , bool calc_sigma
) {
  auto layer = TrkrDefs::getLayer(hitsetkey);
  if (layer <7 || layer > 55) {
    cout << " Error! Trying to compar cluster in layer < 7 or > 55, which is not programmed yet!" << endl;
    return {false, 0.};
  }

  auto clus_T = m_TruthClusterContainer ->findCluster(key_T);
  auto clus_R = m_RecoClusterContainer  ->findCluster(key_R);

  const auto phi_step = m_phistep[layer];

  const auto phi_T    = clus_T->getPosition(0);
  // const auto phisize_T = clus_T->getPhiSize(); // not currently used for matching
  const auto phi_R    = clus_R->getPosition(0);
  const auto phisize_R = clus_R->getPhiSize();

  const auto z_T      = clus_T->getPosition(1);
  // const auto zsize_T   = clus_T->getZSize(); // not currently used for matching
  const auto z_R      = clus_R->getPosition(1);
  const auto zsize_R  = clus_R->getZSize();

  const double phi_delta = abs_dphi (phi_T,phi_R);
  const double z_delta   = fabs (z_T-z_R);

  const double phi_stat = (m_cluster_nphiwidths * phi_step * phisize_R );
  const double z_stat   = (m_cluster_nzwidths   * m_zstep  * zsize_R   );

  if ( phi_delta > phi_stat || z_delta > z_stat ) return { false, 0. };
  else if (!calc_sigma) {
    return { true, 0. };
  } else {
    return { true, phi_delta/phi_stat +z_delta/z_stat };
  }
}

float TruthRecoTrackMatching::sigma_CompMatchClusters(PossibleMatch& match) {
  auto index_true = match[PM_idtrue];
  auto id_reco = match[PM_idreco];

  auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[index_true];
  if (!truth_track) return std::numeric_limits<float>::max();
  std::map<TrkrDefs::hitsetkey,TrkrDefs::cluskey> truth_keys; // truth cluster keys
  for (auto& key : truth_track->getClusters()) truth_keys[TrkrDefs::getHitSetKeyFromClusKey(key)] = key;


  SvtxTrack* reco_track = m_SvtxTrackMap->get(id_reco);
  if (!reco_track) return std::numeric_limits<float>::max();

  auto tpcseed = reco_track->get_tpc_seed();
  if (!tpcseed)    return std::numeric_limits<float>::max();

  double n_matches = 0.; // get the mean match values
  double sum_diff  = 0.;

  for (auto r_key = tpcseed->begin_cluster_keys(); r_key!=tpcseed->end_cluster_keys(); ++r_key) {
    auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(*r_key);
    if (truth_keys.count(hitsetkey) == 0) continue;
    auto comp_val = compare_cluster_pair(truth_keys[hitsetkey], *r_key, hitsetkey, true);

    if (comp_val.first) {
      n_matches += 1.;
      sum_diff += comp_val.second;
    }
  }
  return sum_diff/n_matches;
}

inline bool TruthRecoTrackMatching::skip_match( PossibleMatch& match )
{
  return 
     at_nmax_id_reco    (match[PM_idreco])
  || at_nmax_index_true (match[PM_idreco]);
}

// ----------------------------------------------
// functions for permissible matching pt 
// Currently not used, therefore tracks of any pT
// will be compared against each other. 
// If wanted, will have to be updated with 
// same values in the future.
// ----------------------------------------------
float TruthRecoTrackMatching::delta_outer_pt(float pt) const {
  return -1. + pt*0.; //10. + 0.01*pt;
}
float TruthRecoTrackMatching::delta_inner_pt(float pt) const {
  return -1. + pt*0.; //10. + 0.01*pt;
}

bool TruthRecoTrackMatching::at_nmax_index_true(unsigned short index_true) {
  if (m_nmatched_index_true.find(index_true) == m_nmatched_index_true.end()) {
    return false;
  }
  return m_nmatched_index_true[index_true] >= m_max_nreco_per_truth;
}

bool TruthRecoTrackMatching::at_nmax_id_reco(unsigned short id_reco) {
  if (m_nmatched_id_reco->find(id_reco) == m_nmatched_id_reco->end()) {
    return false;
  }
  return (*m_nmatched_id_reco)[id_reco] >= m_max_ntruth_per_reco;
}

void TruthRecoTrackMatching::add_match_eval(unsigned short id_reco, unsigned short id_true,
    std::map<TrkrDefs::cluskey,TrkrDefs::cluskey> cluskey_map)
{
  if ( !(m_TrackEvalContainer && m_ActsGeometry) )
    return;

  SvtxTrack* track = m_SvtxTrackMap->get(id_reco);
  TrackEvaluationContainerv1::TrackStruct track_struct;

  track_struct.charge = track->get_charge();
  track_struct.nclusters = track->size_cluster_keys();
  track_struct.chisquare = track->get_chisq();
  track_struct.ndf = track->get_ndf();

  track_struct.x = track->get_x();
  track_struct.y = track->get_y();
  track_struct.z = track->get_z();
  track_struct.r = get_r( track_struct.x, track_struct.y );
  track_struct.phi = std::atan2( track_struct.y, track_struct.x );

  track_struct.px = track->get_px();
  track_struct.py = track->get_py();
  track_struct.pz = track->get_pz();
  track_struct.pt = get_pt( track_struct.px, track_struct.py );
  track_struct.p = get_p( track_struct.px, track_struct.py, track_struct.pz );
  track_struct.eta = get_eta( track_struct.p, track_struct.pz );

  // get particle
  auto particle = m_PHG4TruthInfoContainer->GetParticle(id_true);
  if (particle)
  {
    PHG4VtxPoint* vtx  = m_PHG4TruthInfoContainer->GetVtx(particle->get_vtx_id());
    track_struct.embed = m_PHG4TruthInfoContainer->isEmbeded(particle->get_primary_id());
    track_struct.is_primary = particle->get_parent_id() == 0;
    track_struct.pid = particle->get_pid();
    track_struct.gtrackID = particle->get_track_id();
    track_struct.truth_x = vtx->get_x();
    track_struct.truth_y = vtx->get_y();
    track_struct.truth_z = vtx->get_z();
    track_struct.truth_t = vtx->get_t();
    track_struct.truth_px = particle->get_px();
    track_struct.truth_py = particle->get_py();
    track_struct.truth_pz = particle->get_pz();
    track_struct.truth_pt = get_pt( track_struct.truth_px, track_struct.truth_py );
    track_struct.truth_p = get_p( track_struct.truth_px, track_struct.truth_py, track_struct.truth_pz );
    track_struct.truth_eta = get_eta( track_struct.truth_p, track_struct.truth_pz );
  }

  // loop over clusters
  TrackFitUtils::position_vector_t xy_pts;
  for (const auto& key_pair : cluskey_map)
  {
    auto key_R = key_pair.first;
    auto key_T = key_pair.second;
    auto clus_R = m_RecoClusterContainer->findCluster(key_R);
    auto clus_T = m_TruthClusterContainer->findCluster(key_T);
    if ( !(clus_R && clus_T) )
    {
      std::cout << "TruthRecoTrackMatching::add_match_eval - unable to find cluster for key " << clus_R << ":" << clus_T << std::endl;
      continue;
    }

    // create new cluster struct
    Acts::Vector3 global;
    global = m_ActsGeometry->getGlobalPosition(key_R, clus_R);
    TrackEvaluationContainerv1::ClusterStruct cluster_struct;
    cluster_struct.key = key_R;
    cluster_struct.layer = TrkrDefs::getLayer(key_R);
    cluster_struct.x = global.x();
    cluster_struct.y = global.y();
    cluster_struct.z = global.z();
    cluster_struct.r = get_r( cluster_struct.x, cluster_struct.y );
    cluster_struct.phi = std::atan2( cluster_struct.y, cluster_struct.x );

    // truth information
    global = m_ActsGeometry->getGlobalPosition(key_T, clus_T);
    cluster_struct.truth_x = global.x();
    cluster_struct.truth_y = global.y();
    cluster_struct.truth_z = global.z();
    cluster_struct.truth_r = get_r( cluster_struct.truth_x, cluster_struct.truth_y );
    cluster_struct.truth_phi = std::atan2( cluster_struct.truth_y, cluster_struct.truth_x );

    // add to track
    track_struct.clusters.push_back( cluster_struct );

    // add to fitting points
    xy_pts.push_back(std::make_pair(cluster_struct.x, cluster_struct.y));
  }

  // calculate R, X0, and Y0 values
  auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(xy_pts);
  track_struct.R = R;
  track_struct.X0 = X0;
  track_struct.Y0 = Y0;

  m_TrackEvalContainer->addTrack( track_struct );
}
