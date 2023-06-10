#include "TruthRecoTrackMatching.h"

#include "g4evaltools.h"

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>

#include <g4tracking/EmbRecoMatch.h>
#include <g4tracking/EmbRecoMatchContainer.h>
#include <g4tracking/EmbRecoMatchContainerv1.h>
#include <g4tracking/EmbRecoMatchv1.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

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

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>

// To change:
// Make the definition of matching clusters to be if the truth cluster center is withing 1/2 of width of the reco track center

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

TruthRecoTrackMatching::TruthRecoTrackMatching(const unsigned short _nmincluster_match, const float _nmincluster_ratio, const float _cutoff_dphi, const float _same_dphi, const float _cutoff_deta, const float _same_deta, const float _cluster_nzwidths, const float _cluster_nphiwidths, const unsigned short _max_nreco_per_truth, const unsigned short _max_ntruth_per_reco)
  : m_cluster_comp{_cluster_nphiwidths, _cluster_nzwidths}
  , m_nmincluster_match{_nmincluster_match}  // minimum number of clusters to match, default=4
  , m_nmincluster_ratio{_nmincluster_ratio}  // minimum ratio to match a track, default=0.
                                             // -- Track Kinematic Cuts to match --
  , m_cutoff_dphi{_cutoff_dphi}              // maximum |dphi|=|phi_reco-phi_truth| to search for match
  , m_same_dphi{_same_dphi}                  // all tracks in this |dphi| must be tested for matches
  , m_cutoff_deta{_cutoff_deta}              // same as m_cutoff_dphi for deta
  , m_same_deta{_same_deta}                  // same as m_same_dphi for deta
                                             // cluster matching widths (how close the truth center must be reco center)
                                             // number of truth tracks allowed matched per reco track, and v. versa
  , m_max_nreco_per_truth{_max_nreco_per_truth}
  , m_max_ntruth_per_reco{_max_ntruth_per_reco}
{
  m_cluscntr.set_comparer(&m_cluster_comp);
  if (Verbosity() > 50)
  {
    std::cout << " Starting TruthRecoTrackMatching.cc " << std::endl;
  }
}

int TruthRecoTrackMatching::InitRun(PHCompositeNode* topNode)  //`
{
  if (Verbosity() > 10)
  {
    topNode->print();
  }
  auto init_status = m_cluster_comp.init(topNode);
  if (init_status == Fun4AllReturnCodes::ABORTRUN)
  {
    return init_status;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_nmatched_id_reco = &(m_EmbRecoMatchContainer->map_nTruthPerReco());
  m_nmatched_id_true = &(m_EmbRecoMatchContainer->map_nRecoPerTruth());
  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthRecoTrackMatching::process_event(PHCompositeNode* topnode)
{
  if (topnode == nullptr)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (Verbosity() > 1000)
  {
    topnode->print();  // perhaps not needed
  }

  // cleanup output container
  if (m_TrackEvalContainer) {
    m_TrackEvalContainer->Reset();
    m_TrackEvalContainer->clearEvents();
    m_TrackEvalContainer->clearClusters();
    m_TrackEvalContainer->clearTracks();
  }
  m_nmatched_index_true.clear();

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
  if (Verbosity() > 60)
  {
    std::cout << "reco tracks size: " << m_SvtxTrackMap->size() << std::endl;
  }

  recoData.clear();

  for (auto& reco : *m_SvtxTrackMap)
  {
    auto index_reco = reco.first;
    auto track = reco.second;
    recoData.push_back({track->get_phi(), track->get_eta(), track->get_pt(), index_reco});
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
  auto wrap_from_start = std::upper_bound(recoData.begin(),
                                          recoData.end(), (-M_PI + m_cutoff_dphi), CompRECOtoPhi());
  for (auto iter = recoData.begin(); iter != wrap_from_start; ++iter)
  {
    auto entry = *iter;                                              // make a new copy to wrap to the other side of recoData
    std::get<RECOphi>(entry) = std::get<RECOphi>(entry) + 2 * M_PI;  // put the new copy on the other end
    wraps.push_back(entry);
  }
  for (auto E : wraps)
  {
    recoData.push_back(E);
  }

  if (Verbosity() > 70)
  {
    std::cout << " ****************************************** " << std::endl;
    std::cout << " This is the RECO map of tracks to match to " << std::endl;
    std::cout << " ****************************************** " << std::endl;
    for (auto& E : recoData)
    {
      std::cout << Form(" id:%2i  (phi:eta:pt) (%5.2f:%5.2f:%5.2f)", std::get<RECOid>(E),
                        std::get<RECOphi>(E), std::get<RECOeta>(E), std::get<RECOpt>(E))
                << std::endl;
    }
    std::cout << " ****end*listing*********************** " << std::endl;
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

  std::vector<std::pair<unsigned short, unsigned short>> outerbox_pairs{};
  std::vector<std::pair<unsigned short, unsigned short>> innerbox_pairs{};

  if (topnode == nullptr)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (Verbosity() > 70)
  {
    std::cout << "Number of truth tracks: " << m_TrkrTruthTrackContainer->getMap().size() << std::endl;
    for (auto& _pair : m_TrkrTruthTrackContainer->getMap())
    {
      auto& track = _pair.second;
      std::cout << Form(" id:%2i  (phi:eta:pt) (%5.2f:%5.2f:%5.2f nclus: %i)", track->getTrackid(),
                        track->getPhi(), track->getPseudoRapidity(), track->getPt(),
                        (int) track->getClusters().size())
                << std::endl;
    }
    std::cout << " ----end-listing-truth-tracks---------- " << std::endl;
  }

  for (auto _pair : m_TrkrTruthTrackContainer->getMap())
  {
    auto id_true = _pair.first;
    auto track = _pair.second;

    auto match_indices = find_box_matches(track->getPhi(),
                                          track->getPseudoRapidity(), track->getPt());

    // keep track of all truth tracks (by track-id) which has been matched
    for (auto& id_reco : match_indices.first)
    {
      innerbox_pairs.emplace_back(id_true, id_reco);
    }
    for (auto& id_reco : match_indices.second)
    {
      outerbox_pairs.emplace_back(id_true, id_reco);
    }

    if (Verbosity() > 80)
    {
      std::cout << Form(" T(%i)  find box(%5.2f:%5.2f:%5.2f)",
                        (int) track->getTrackid(), track->getPhi(),
                        track->getPseudoRapidity(), track->getPt());
      for (auto& id_reco : match_indices.first)
      {
        std::cout << "->IB(" << id_reco << ")";
      }
      for (auto& id_reco : match_indices.second)
      {
        std::cout << "->OB(" << id_reco << ")";
      }
      std::cout << std::endl;
    }
  }

  if (Verbosity() > 100)
  {
    std::cout << "Innerbox pairs" << std::endl;
  }
  match_tracks_in_box(innerbox_pairs);
  if (Verbosity() > 100)
  {
    std::cout << "Outerbox pairs" << std::endl;
  }
  match_tracks_in_box(outerbox_pairs);

  // fill the list of all truth track ids that are not matched
  for (auto _pair : m_TrkrTruthTrackContainer->getMap())
  {
    auto id_true = _pair.first;
    m_EmbRecoMatchContainer->checkfill_idsTruthUnmatched(id_true);
  }

  if (Verbosity() > 50)
  {
    std::cout << " --START--print-all-stored-matches--" << std::endl;
    std::cout << " --0-- Printing all matches stored (start)" << std::endl;
    // re-print all the tracks with the matches with the fit values
    for (auto match : m_EmbRecoMatchContainer->getMatches())
    {
      std::cout << Form(" Match id(%2i->%2i) nClusMatch-nClusTrue-nClusReco (%2i:%2i:%2i)",
                        match->idTruthTrack(), match->idRecoTrack(),
                        match->nClustersMatched(), match->nClustersTruth(), match->nClustersReco())
                << std::endl;
    }
    std::cout << " --END--print-all-stored-matches----" << std::endl;
  }

  if (m_write_diag)
  {
    fill_tree();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthRecoTrackMatching::End(PHCompositeNode* /*topNode*/)
{
  TFile* s_current = gDirectory->GetFile();
  if (m_write_diag)
  {
    m_diag_file->cd();
    G4Eval::write_StringToTFile(
        "trk_match_sel",
        Form("  Matching criteria:\n"
             "  min. clusters to match:  %i\n"
             "  min. clust. match ratio: %4.2f"
             "  dphi small window:       %4.2f"
             "  dphi large windows:      %4.2f"
             "  deta small window:       %4.2f"
             "  deta large window:       %4.2f"
             "  nmax phg4 matches per svtx:  %i"
             "  nmax svtx matches per phg4:  %i",
             m_nmincluster_match, m_nmincluster_ratio, m_same_dphi, m_cutoff_dphi, m_same_deta, m_cutoff_deta, m_max_ntruth_per_reco, m_max_nreco_per_truth));
    m_diag_tree->Write();
    m_diag_file->Save();
    m_diag_file->Close();
    if (s_current != nullptr)
    {
      s_current->cd();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//--------------------------------------------------
// Internal functions
//--------------------------------------------------

int TruthRecoTrackMatching::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // Initiailize the modules
  m_PHG4TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_PHG4TruthInfoContainer)
  {
    std::cout << "Could not locate G4TruthInfo node when running "
              << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_SvtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_SvtxTrackMap)
  {
    std::cout << "Could not locate SvtxTrackMap node when running "
              << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_TruthClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode,
                                                                     "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_TruthClusterContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHCLUSTERCONTAINER node when "
              << "running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_RecoClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_RecoClusterContainer)
  {
    std::cout << "Could not locate TRKR_CLUSTER node when running "
              << "\"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_TrkrTruthTrackContainer = findNode::getClass<TrkrTruthTrackContainer>(topNode,
                                                                          "TRKR_TRUTHTRACKCONTAINER");
  if (!m_TrkrTruthTrackContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHTRACKCONTAINER node when "
              << "running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_PHG4TpcCylinderGeomContainer =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_PHG4TpcCylinderGeomContainer)
  {
    std::cout << "Could not locate CYLINDERCELLGEOM_SVTX node when "
              << "running \"TruthRecoTrackMatching\" module." << std::endl;
    /* std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl; */
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // note that layers 0-6, and > 55, don't work
  /* for (int layer=7; layer<55; ++layer) { */
  /*   PHG4TpcCylinderGeom *layergeom = m_PHG4TpcCylinderGeomContainer->GetLayerCellGeom(layer); */
  /*   if (layer==7) m_zstep = layergeom->get_zstep(); */
  /*   m_phistep[layer] = layergeom->get_phistep(); */
  /* } */
  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode,
                                                                      "TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer)
  {
    PHNodeIterator dstiter(dstNode);

    auto DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    m_EmbRecoMatchContainer = new EmbRecoMatchContainerv1;
    auto newNode = new PHIODataNode<PHObject>(m_EmbRecoMatchContainer,
                                              "TRKR_EMBRECOMATCHCONTAINER", "PHObject");
    DetNode->addNode(newNode);

    m_TrackEvalContainer = new TrackEvaluationContainerv1;
    auto newNode2 = new PHIODataNode<PHObject>(m_TrackEvalContainer, "TrackMatchEvalContainer", "PHObject");
    DetNode->addNode(newNode2);
  }

  m_ActsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  return Fun4AllReturnCodes::EVENT_OK;
}

std::pair<std::vector<unsigned short>, std::vector<unsigned short>>
TruthRecoTrackMatching::find_box_matches(float truth_phi, float truth_eta, float truth_pt)
{
  // sort through the recoData to find:
  //     inner_box : possible matches withing m_same_d{phi,eta,pt}
  //     outer_box : (assumed to be strictly larger than, and therefore contain, inner_box, in m_cutoff_d{phi, eta,pt}
  // recoData is already sorted by eta, but as the boxes "zoom in" on acceptable matches, it gets
  // sorted by eta and (if applicable) pT as well. `mix_pair` keeps track of the range of recoData that has
  // been sorted so that it can be resorted to phi order.
  RECO_pair_iter mix_pair{recoData.begin(), recoData.end()};

  mix_pair.first = std::lower_bound(mix_pair.first, mix_pair.second,
                                    truth_phi - m_cutoff_dphi, CompRECOtoPhi());

  mix_pair.second = std::upper_bound(mix_pair.first, mix_pair.second,
                                     truth_phi + m_cutoff_dphi, CompRECOtoPhi());

  if (mix_pair.first == mix_pair.second)
  {
    // there are no possible phi_matches; return nothing
    return {{}, {}};
  }

  // There are tracks within the outer_box phi range;
  // now re-shuffle the phi outerbox range for eta and see if there
  // are tracks in outer_box-eta range
  std::sort(mix_pair.first, mix_pair.second, CompRECOtoEta());

  RECO_pair_iter outer_box = mix_pair;
  outer_box.first = std::lower_bound(mix_pair.first, mix_pair.second,
                                     truth_eta - m_cutoff_deta, CompRECOtoEta());

  outer_box.second = std::lower_bound(outer_box.first, outer_box.second,
                                      truth_eta + m_cutoff_deta, CompRECOtoEta());

  if (outer_box.first == outer_box.second)
  {
    // there are no eta_matches in outerbox; resort mix_pair and return nothing
    std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());
    return {{}, {}};
  }

  // if pt is specified, restrict pt_box to the given range
  RECO_pair_iter inner_box = outer_box;
  const float _delta_outer_pt = delta_outer_pt(truth_pt);
  const float _delta_inner_pt = delta_inner_pt(truth_pt);
  if (_delta_outer_pt > 0)
  {
    std::sort(outer_box.first, outer_box.second, CompRECOtoPt());
    outer_box.first = std::lower_bound(outer_box.first, outer_box.second, truth_pt - _delta_outer_pt, CompRECOtoPt());
    outer_box.second = std::upper_bound(outer_box.first, outer_box.second, truth_pt + _delta_outer_pt, CompRECOtoPt());

    // if not outer_box, resort and return nothing
    if (outer_box.first == outer_box.second)
    {
      // there are no eta_matches in outerbox; resort mix_pair and return nothing
      std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());
      return {{}, {}};
    }

    // now calculate the inner box -- this time sorting in reverse -- pT, then eta, then phi
    inner_box = outer_box;
    if (_delta_inner_pt > 0)
    {
      // start the inner pair -- the outerbox is already sorted by pT
      inner_box.first = std::lower_bound(inner_box.first,
                                         inner_box.second, truth_pt - _delta_inner_pt, CompRECOtoPt());
      inner_box.second = std::upper_bound(inner_box.first,
                                          inner_box.second, truth_pt + _delta_inner_pt, CompRECOtoPt());
    }

    // go back to sorted by eta
    std::sort(inner_box.first, inner_box.second, CompRECOtoEta());
  }

  // At this point we know that we have outer_box matches
  // Finish checking if there are any possible inner_box matches

  // check for inner_box eta matches
  if (inner_box.first != inner_box.second)
  {
    inner_box.first = std::lower_bound(inner_box.first, inner_box.second,
                                       truth_eta - m_same_deta, CompRECOtoEta());

    inner_box.second = std::lower_bound(inner_box.first, inner_box.second,
                                        truth_eta + m_same_deta, CompRECOtoEta());
  }

  // check for inner_box phi matches
  if (inner_box.first != inner_box.second)
  {
    std::sort(inner_box.first, inner_box.second, CompRECOtoPhi());

    inner_box.first = std::lower_bound(inner_box.first, inner_box.second,
                                       truth_phi - m_same_dphi, CompRECOtoPhi());

    inner_box.second = std::lower_bound(inner_box.first, inner_box.second,
                                        truth_phi + m_same_dphi, CompRECOtoPhi());
  }

  // At this point there are definitely outer_box matches, and
  // possible inner_box matches.
  // Return all these possible matches, evaluating all inner_box_matches
  std::vector<unsigned short> inner_box_matches, outer_box_matches;

  // fill inner_box_matches
  for (auto iter = inner_box.first; iter != inner_box.second; ++iter)
  {
    inner_box_matches.push_back(std::get<RECOid>(*iter));
  }
  // fill inner_box_matches
  for (auto iter = outer_box.first; iter != inner_box.first; ++iter)
  {
    outer_box_matches.push_back(std::get<RECOid>(*iter));
  }
  for (auto iter = inner_box.second; iter != outer_box.second; ++iter)
  {
    outer_box_matches.push_back(std::get<RECOid>(*iter));
  }

  // resort recoData back to phi order
  std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());

  // return the box matches
  return std::make_pair(inner_box_matches, outer_box_matches);
}

void TruthRecoTrackMatching::match_tracks_in_box(
    std::vector<std::pair<unsigned short, unsigned short>>& box_pairs  // possible matches by proximity in eta, phi, pT
)
{
  if (box_pairs.size() == 0)
  {
    return;
  }

  std::sort(box_pairs.begin(), box_pairs.end());  // sorted by first index_true, then id_reco
  std::vector<PossibleMatch> poss_matches;
  vector<std::map<TrkrDefs::cluskey,TrkrDefs::cluskey>> v_clusmap; // key_reco, key_true
  unsigned short index_clusmap = 0;

  auto ipair = box_pairs.begin();  // save current examined pair for sorting/unsorting purposes
  while (ipair != box_pairs.end())
  {
    auto id_true = ipair->first;

    // See if this track already has the maximum number of matches. If it does, then skip all pairs using this truth-track-index;
    if (at_nmax_index_true(id_true))
    {
      while (ipair != box_pairs.end())
      {
        ++ipair;
        if (ipair->first != id_true)
        {
          break;
        }
      }
      continue;
    }

    // make all possible reco matches matched for this track
    /* std::map<TrkrDefs::hitsetkey,TrkrDefs::cluskey> truth_keys; */
    m_cluscntr.addClusKeys(m_TrkrTruthTrackContainer->getTruthTrack(id_true));
    // add the truth keys into the track counter
    /* auto truth_track = m_TrkrTruthTrackContainer->getTruthTrack(id_true); */
    /* for (auto& key : truth_track->getClusters()) { */
    /*   truth_keys[TrkrDefs::getHitSetKeyFromClusKey(key)] = key; */
    /* } */
    /* unsigned short nclus_true = truth_keys.size(); */

    while (ipair != box_pairs.end())
    {  // Loop over all possible matches only for this index_true
      //(subsequent true tracks will be caught on following loops)
      if (ipair->first != id_true)
      {
        break;  // never true on first iteration
      }
      if (at_nmax_id_reco(ipair->second))
      {
        ++ipair;
        continue;
      }  // make sure the reco-track isn't alread matched
      // ok, make a possible match: compare the clusters in the truth track and the reco track
      SvtxTrack* reco_track = m_SvtxTrackMap->get(ipair->second);
      m_cluscntr.addClusKeys(reco_track);
      m_cluscntr.find_matches();

      /* unsigned short nclus_match   = 0; // fill in the comparison loop */
      /* unsigned short nclus_nomatch = 0; // number of reco and truth cluster
       *             // that share a hitsetkey, but still fair matching criteria */
      /* unsigned short nclus_reco    = 0; // count in the comparison loop */
      // do the comparison of the tracks
      /* for (auto reco_ckey : G4Eval::ClusKeyIter(reco_track)) { */
      /*   ++nclus_reco; */
      /*   auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(reco_ckey); */
      /*   if (truth_keys.count(hitsetkey) != 0) { */
      /*     // reco and truth cluster are in same hitsetkey-indexed subsurface. */
      /*     // See if they match (++nclus_match) or not (++nclus_nomatch) */
      /*     if (m_cluster_comp(truth_keys[hitsetkey], reco_ckey).first) { */
      /*       ++nclus_match; */
      /*     } else { */
      /*       ++nclus_nomatch; */
      /*     } */
      /*   } */
      /* } */
      unsigned short nclus_match = m_cluscntr.phg4_n_matched();
      unsigned short nclus_true = m_cluscntr.phg4_nclus();
      unsigned short nclus_reco = m_cluscntr.svtx_nclus();
      unsigned short nclus_nomatch = (int) (m_cluscntr.svtx_keys.size() - nclus_match);

      if (Verbosity() > 100)
      {
        auto truth_track = m_TrkrTruthTrackContainer->getTruthTrack(id_true);
        std::cout << Form(
                         "possmatch:(phi,eta,pT:id) true(%5.2f,%5.2f,%4.2f:%2i) reco(%5.2f,%5.2f,%4.2f:%2i) "
                         "nCl(match:true:reco:nomatch)(%2i-%2i-%2i-%2i)",
                         truth_track->getPhi(), truth_track->getPseudoRapidity(), truth_track->getPt(),
                         (int) truth_track->getTrackid(),
                         reco_track->get_phi(), reco_track->get_eta(), reco_track->get_pt(),
                         ipair->second,
                         (int) nclus_match, (int) nclus_true, (int) nclus_reco, (int) nclus_nomatch)
                  << std::endl;
      }
      if (nclus_match >= m_nmincluster_match && (static_cast<float>(nclus_match) / nclus_true >= m_nmincluster_ratio))
      {
        poss_matches.push_back(
            {nclus_match, nclus_true, nclus_reco,
             ipair->first, ipair->second});
        v_clusmap.push_back(m_cluscntr.svtx_phg4_keys);
        index_clusmap++;
      }
      ++ipair;
    }
  }

  // add all possible matches started for the largest PM_nmatch (the top)
  //  for groups of ties of PM_nmatch, sort by PM_ntrue (from smallest)
  //    for groups of ties (PM_nmatch, PM_ntrue), go by smallest PM_nreco
  //       for groups of ties (PM_nmatch, PM_ntrue, PM_nreco), do a detailed sum diff in the deltas
  std::sort(poss_matches.begin(), poss_matches.end(), SortPossibleMatch());

  if (Verbosity() > 200)
  {
    std::cout << " All chosen possible matches (" << poss_matches.size() << ") track pairs  (nClMatched-nClTrue-nClReco : idTrue-idReco)" << std::endl;
    int i{0};
    for (auto match : poss_matches)
    {
      /* auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[std::get<PM_idtrue>(match)]; */
      /* auto index_trut = truth_track->getTrackid(); */
      std::cout << Form(" pair(%2i):  %2i-%2i-%2i-<%2i>-%2i ", i++, (int) std::get<PM_nmatch>(match), (int) std::get<PM_ntrue>(match), (int) std::get<PM_nreco>(match), (int) std::get<PM_idtrue>(match), (int) std::get<PM_idreco>(match)) << std::endl;
    }
  }

  /* std::set<int> matched_idreco, matched_idtrue; */
  auto iter = poss_matches.begin();
  while (iter != poss_matches.end())
  {
    if (skip_match(*iter))
    {
      ++iter;
      continue;
    }
    std::vector<std::pair<float, int>> sigma_metric = {{0., 0}};
    int n_sigma = 0;
    auto iter_same = iter + 1;  // iterator to look forward from first point and get detailed comparisons for all possibly matched tracks
    while (
        iter_same != poss_matches.end() && (*iter_same)[PM_nmatch] == (*iter)[PM_nmatch] && (*iter_same)[PM_ntrue] == (*iter)[PM_ntrue] && (*iter_same)[PM_nreco] == (*iter)[PM_nreco])
    {
      ++n_sigma;
      if (n_sigma == 1)
      {
        sigma_metric[0].first = sigma_CompMatchClusters(*iter);
      }
      sigma_metric.emplace_back(sigma_CompMatchClusters(*iter_same), n_sigma);
      ++iter_same;
    }
    std::sort(sigma_metric.begin(), sigma_metric.end());

    bool first = true;
    for (auto& sigma : sigma_metric)
    {
      if (first)
      {
        first = false;
      }
      else if (skip_match(*(iter + sigma.second)))
      {
        continue;
      }
      auto match = *(iter + sigma.second);
      auto id_true = match[PM_idtrue];
      auto id_reco = match[PM_idreco];
      auto save_match = new EmbRecoMatchv1(id_true, id_reco,
                                           match[PM_ntrue], match[PM_nreco], match[PM_nmatch]);
      m_EmbRecoMatchContainer->addMatch(save_match);
      add_match_eval(id_reco, id_true, v_clusmap[match[PM_clusmap]]);

      if (m_nmatched_index_true.find(id_true) == m_nmatched_index_true.end())
      {
        m_nmatched_index_true[id_true] = 1;
      }
      else
      {
        m_nmatched_index_true[id_true] += 1;
      }
    }
    iter += sigma_metric.size();
  }
}

// ----------------------------------------
// convenience function
// ----------------------------------------
inline float TruthRecoTrackMatching::abs_dphi(float phi0, float phi1)
{
  float dphi = std::fabs(phi0 - phi1);
  while (dphi > M_PI)
  {
    dphi = std::fabs(dphi - 2 * M_PI);
  }
  return dphi;
}

float TruthRecoTrackMatching::sigma_CompMatchClusters(PossibleMatch& match)
{
  auto id_true = match[PM_idtrue];
  auto id_reco = match[PM_idreco];

  m_cluscntr.addClusKeys(m_TrkrTruthTrackContainer->getTruthTrack(id_true));
  m_cluscntr.addClusKeys(m_SvtxTrackMap->get(id_reco));

  m_cluscntr.find_matches();

  int n_matched = m_cluscntr.phg4_n_matched();

  if (n_matched)
  {
    return std::numeric_limits<float>::max();
  }
  else
  {
    return m_cluscntr.match_stat / n_matched;
  }
}
/* auto truth_track = m_TrkrTruthTrackContainer->getTruthTrack(id_true); */
/* if (!truth_track) return std::numeric_limits<float>::max(); */
/* std::map<TrkrDefs::hitsetkey,TrkrDefs::cluskey> truth_keys; // truth cluster keys */
/* for (auto& key : truth_track->getClusters()) */
/*   truth_keys[TrkrDefs::getHitSetKeyFromClusKey(key)] = key; */

/* SvtxTrack* reco_track = m_SvtxTrackMap->get(id_reco); */
/* if (!reco_track) return std::numeric_limits<float>::max(); */

/* auto tpcseed = reco_track->get_tpc_seed(); */
/* if (!tpcseed)    return std::numeric_limits<float>::max(); */

/* double n_matches = 0.; // get the mean match values */
/* double sum_diff  = 0.; */

/* for (auto reco_ckey : G4Eval::ClusKeyIter(reco_track)) { */
/* auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(reco_ckey); */
/* if (truth_keys.count(hitsetkey) == 0) continue; */
/* auto comp_val = m_cluster_comp(truth_keys[hitsetkey], reco_ckey); */

/* if (comp_val.first) { */
/* n_matches += 1.; */
/* sum_diff += comp_val.second; */
/* } */
/* } */
/* return sum_diff/n_matches; */
/* } */

inline bool TruthRecoTrackMatching::skip_match(PossibleMatch& match)
{
  return at_nmax_id_reco(match[PM_idreco]) || at_nmax_index_true(match[PM_idreco]);
}

// ----------------------------------------------
// functions for permissible matching pt
// Currently not used, therefore tracks of any pT
// will be compared against each other.
// If wanted, will have to be updated with
// same values in the future.
// ----------------------------------------------
float TruthRecoTrackMatching::delta_outer_pt(float pt) const
{
  return -1. + pt * 0.;  // 10. + 0.01*pt;
}
float TruthRecoTrackMatching::delta_inner_pt(float pt) const
{
  return -1. + pt * 0.;  // 10. + 0.01*pt;
}

bool TruthRecoTrackMatching::at_nmax_index_true(unsigned short index_true)
{
  if (m_nmatched_index_true.find(index_true) == m_nmatched_index_true.end())
  {
    return false;
  }
  return m_nmatched_index_true[index_true] >= m_max_nreco_per_truth;
}

bool TruthRecoTrackMatching::at_nmax_id_reco(unsigned short id_reco)
{
  if (m_nmatched_id_reco->find(id_reco) == m_nmatched_id_reco->end())
  {
    return false;
  }
  return (*m_nmatched_id_reco)[id_reco] >= m_max_ntruth_per_reco;
}

void TruthRecoTrackMatching::set_diagnostic_file(const std::string& file_name)
{
  m_write_diag = true;
  TFile* s_current = gDirectory->GetFile();
  m_diag_file = new TFile(file_name.c_str(), "recreate");
  // write out the cuts on this set of data
  m_diag_file->cd();
  m_diag_tree = new TTree("T", "Tree of Reco and True Clusters");

  m_diag_tree->Branch("event", &m_event);

  m_diag_tree->Branch("trkid_reco_matched", &m_trkid_reco_matched);
  m_diag_tree->Branch("itrk_reco_matched", &m_cnt_reco_matched);
  m_diag_tree->Branch("i0_reco_matched", &m_i0_reco_matched);
  m_diag_tree->Branch("i1_reco_matched", &m_i1_reco_matched);
  m_diag_tree->Branch("layer_reco_matched", &m_layer_reco_matched);
  m_diag_tree->Branch("x_reco_matched", &m_x_reco_matched);
  m_diag_tree->Branch("y_reco_matched", &m_y_reco_matched);
  m_diag_tree->Branch("z_reco_matched", &m_z_reco_matched);

  m_diag_tree->Branch("trkid_reco_notmatched", &m_trkid_reco_notmatched);
  m_diag_tree->Branch("itrk_reco_notmatched", &m_cnt_reco_notmatched);
  m_diag_tree->Branch("i0_reco_notmatched", &m_i0_reco_notmatched);
  m_diag_tree->Branch("i1_reco_notmatched", &m_i1_reco_notmatched);
  m_diag_tree->Branch("layer_reco_notmatched", &m_layer_reco_notmatched);
  m_diag_tree->Branch("x_reco_notmatched", &m_x_reco_notmatched);
  m_diag_tree->Branch("y_reco_notmatched", &m_y_reco_notmatched);
  m_diag_tree->Branch("z_reco_notmatched", &m_z_reco_notmatched);

  m_diag_tree->Branch("trkid_true_matched", &m_trkid_true_matched);
  m_diag_tree->Branch("itrk_true_matched", &m_cnt_true_matched);
  m_diag_tree->Branch("i0_true_matched", &m_i0_true_matched);
  m_diag_tree->Branch("i1_true_matched", &m_i1_true_matched);
  m_diag_tree->Branch("layer_true_matched", &m_layer_true_matched);
  m_diag_tree->Branch("x_true_matched", &m_x_true_matched);
  m_diag_tree->Branch("y_true_matched", &m_y_true_matched);
  m_diag_tree->Branch("z_true_matched", &m_z_true_matched);

  m_diag_tree->Branch("trkid_true_notmatched", &m_trkid_true_notmatched);
  m_diag_tree->Branch("itrk_true_notmatched", &m_cnt_true_notmatched);
  m_diag_tree->Branch("i0_true_notmatched", &m_i0_true_notmatched);
  m_diag_tree->Branch("i1_true_notmatched", &m_i1_true_notmatched);
  m_diag_tree->Branch("layer_true_notmatched", &m_layer_true_notmatched);
  m_diag_tree->Branch("x_true_notmatched", &m_x_true_notmatched);
  m_diag_tree->Branch("y_true_notmatched", &m_y_true_notmatched);
  m_diag_tree->Branch("z_true_notmatched", &m_z_true_notmatched);

  if (s_current != nullptr)
  {
    s_current->cd();
  }
}

void TruthRecoTrackMatching::clear_branch_vectors()
{
  m_trkid_reco_matched.clear();
  m_i0_reco_matched.clear();
  m_i1_reco_matched.clear();
  m_layer_reco_matched.clear();
  m_x_reco_matched.clear();
  m_y_reco_matched.clear();
  m_z_reco_matched.clear();

  m_trkid_reco_notmatched.clear();
  m_i0_reco_notmatched.clear();
  m_i1_reco_notmatched.clear();
  m_layer_reco_notmatched.clear();
  m_x_reco_notmatched.clear();
  m_y_reco_notmatched.clear();
  m_z_reco_notmatched.clear();

  m_trkid_true_matched.clear();
  m_i0_true_matched.clear();
  m_i1_true_matched.clear();
  m_layer_true_matched.clear();
  m_x_true_matched.clear();
  m_y_true_matched.clear();
  m_z_true_matched.clear();

  m_trkid_true_notmatched.clear();
  m_i0_true_notmatched.clear();
  m_i1_true_notmatched.clear();
  m_layer_true_notmatched.clear();
  m_x_true_notmatched.clear();
  m_y_true_notmatched.clear();
  m_z_true_notmatched.clear();
}

void TruthRecoTrackMatching::fill_tree()
{
  // fill clusters or un-matched truth tracks
  int cnt = 0;
  int itrk = 0;
  for (auto& trkid : m_EmbRecoMatchContainer->ids_TruthUnmatched())
  {
    m_trkid_true_notmatched.push_back(trkid);
    m_i0_true_notmatched.push_back(cnt);
    auto track = m_TrkrTruthTrackContainer->getTruthTrack(trkid);
    for (auto& ckey : track->getClusters())
    {
      auto cluster = m_TruthClusterContainer->findCluster(ckey);
      m_cnt_true_notmatched.push_back(itrk);
      Eigen::Vector3d gloc = m_ActsGeometry->getGlobalPosition(ckey, cluster);
      m_layer_true_notmatched.push_back(TrkrDefs::getLayer(ckey));
      m_x_true_notmatched.push_back(gloc[0]);
      m_y_true_notmatched.push_back(gloc[1]);
      m_z_true_notmatched.push_back(gloc[2]);
      ++cnt;
    }
    m_i1_true_notmatched.push_back(cnt);
    ++itrk;
  }

  // fill clusters of matched truth tracks
  cnt = 0;
  itrk = 0;
  for (auto& trkid : m_EmbRecoMatchContainer->ids_TruthMatched())
  {
    m_trkid_true_matched.push_back(trkid);
    m_i0_true_matched.push_back(cnt);
    auto track = m_TrkrTruthTrackContainer->getTruthTrack(trkid);
    for (auto& ckey : track->getClusters())
    {
      auto cluster = m_TruthClusterContainer->findCluster(ckey);
      m_cnt_true_matched.push_back(itrk);
      Eigen::Vector3d gloc = m_ActsGeometry->getGlobalPosition(ckey, cluster);
      m_layer_true_matched.push_back(TrkrDefs::getLayer(ckey));
      m_x_true_matched.push_back(gloc[0]);
      m_y_true_matched.push_back(gloc[1]);
      m_z_true_matched.push_back(gloc[2]);
      ++cnt;
    }
    m_i1_true_matched.push_back(cnt);
    ++itrk;
  }

  // fill clusters of matched reco tracks
  std::set<unsigned int> set_reco_matched;
  cnt = 0;
  itrk = 0;
  for (auto& trkid : m_EmbRecoMatchContainer->ids_RecoMatched())
  {
    set_reco_matched.insert(trkid);
    m_trkid_reco_matched.push_back(trkid);
    SvtxTrack* reco_track = m_SvtxTrackMap->get(trkid);
    m_i0_reco_matched.push_back(cnt);

    for (auto reco_ckey : G4Eval::ClusKeyIter(reco_track))
    {
      auto cluster = m_RecoClusterContainer->findCluster(reco_ckey);
      Eigen::Vector3d gloc = m_ActsGeometry->getGlobalPosition(reco_ckey, cluster);
      m_layer_reco_matched.push_back(TrkrDefs::getLayer(reco_ckey));
      m_cnt_reco_matched.push_back(itrk);
      m_x_reco_matched.push_back(gloc[0]);
      m_y_reco_matched.push_back(gloc[1]);
      m_z_reco_matched.push_back(gloc[2]);
      ++cnt;
    }
    m_i1_reco_matched.push_back(cnt);
    ++itrk;
  }

  // fill clusters of not matched reco tracks
  cnt = 0;
  itrk = 0;
  for (auto& reco : *m_SvtxTrackMap)
  {
    auto trkid = reco.first;
    if (set_reco_matched.count(trkid))
    {
      continue;
    }
    m_trkid_reco_notmatched.push_back(trkid);
    SvtxTrack* reco_track = m_SvtxTrackMap->get(trkid);
    m_i0_reco_notmatched.push_back(cnt);
    for (auto reco_ckey : G4Eval::ClusKeyIter(reco_track))
    {
      auto cluster = m_RecoClusterContainer->findCluster(reco_ckey);
      Eigen::Vector3d gloc = m_ActsGeometry->getGlobalPosition(reco_ckey, cluster);
      m_layer_reco_notmatched.push_back(TrkrDefs::getLayer(reco_ckey));
      m_cnt_reco_notmatched.push_back(itrk);
      m_x_reco_notmatched.push_back(gloc[0]);
      m_y_reco_notmatched.push_back(gloc[1]);
      m_z_reco_notmatched.push_back(gloc[2]);
      ++cnt;
    }
    m_i1_reco_notmatched.push_back(cnt);
    ++itrk;
  }
  ++m_event;
  m_diag_tree->Fill();
  clear_branch_vectors();
  return;
}

void TruthRecoTrackMatching::add_match_eval(unsigned short id_reco, unsigned short id_true,
    std::map<TrkrDefs::cluskey,TrkrDefs::cluskey> cluskey_map)
{
  if ( !(m_TrackEvalContainer && m_ActsGeometry) )
    return;

  SvtxTrack* track = m_SvtxTrackMap->get(id_reco);
  TrackEvaluationContainerv1::TrackStruct track_struct;

  track_struct.trackID = track->get_id();
  track_struct.charge = track->get_charge();
  track_struct.nclusters = track->size_cluster_keys();
  track_struct.chisquare = track->get_chisq();
  track_struct.ndf = track->get_ndf();

  track_struct.x = track->get_x();
  track_struct.y = track->get_y();
  track_struct.z = track->get_z();
  track_struct.r = get_r( track_struct.x, track_struct.y );

  track_struct.px = track->get_px();
  track_struct.py = track->get_py();
  track_struct.pz = track->get_pz();
  track_struct.pt = get_pt( track_struct.px, track_struct.py );
  track_struct.p = get_p( track_struct.px, track_struct.py, track_struct.pz );
  track_struct.eta = get_eta( track_struct.p, track_struct.pz );
  track_struct.phi = std::atan2( track_struct.py, track_struct.px );

  // get particle
  auto particle = m_PHG4TruthInfoContainer->GetParticle(id_true);
  if (particle)
  {
    track_struct.embed = m_PHG4TruthInfoContainer->isEmbeded(particle->get_primary_id());
    if (_scan_for_embedded && !track_struct.embed) return;
    track_struct.is_primary = particle->get_parent_id() == 0;
    if (_scan_for_primaries && !track_struct.is_primary) return;
    track_struct.pid = particle->get_pid();
    track_struct.gtrackID = particle->get_track_id();
    track_struct.parentID = particle->get_parent_id();

    PHG4VtxPoint* vtx  = m_PHG4TruthInfoContainer->GetVtx(particle->get_vtx_id());
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
    track_struct.truth_phi = std::atan2( track_struct.truth_py, track_struct.truth_px );
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

  for(int li = 0; li < track_struct.m_nCaloLayers; li++)
  {
    SvtxTrack::CAL_LAYER layer = static_cast<SvtxTrack::CAL_LAYER>(li + 1);
    track_struct.cal_energy_3x3[li] = track->get_cal_energy_3x3(layer);
    track_struct.cal_energy_5x5[li] = track->get_cal_energy_5x5(layer);
    track_struct.cal_cluster_e[li] = track->get_cal_cluster_e(layer);
  }

  m_TrackEvalContainer->addTrack( track_struct );
}
