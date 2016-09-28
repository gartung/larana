#ifndef BACKTRACKMATCHERALG_H
#define BACKTRACKMATCHERALG_H
/*!
 * Title:   BackTrack Matcher Algorithim Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for matching MCParticles to reconstructed hits, clusters, and tracks
 * Input:       simb::MCParticle and (vector<recob::Hit> or vector<recob::Cluster> or vector<recob::Track>)
 * Output:      vector<recob::Hit> or vector<recob::Cluster> or recob::Track and matching info
*/

//Framework includes
#include "fhiclcpp/ParameterSet.h"

//LArSoft includes
#include "larsim/MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"

//c++ includes
#include <vector>
#include <utility>

namespace mctrue {
  class BackTrackMatcherAlg;
}

/// Matches Reco objects to MCParticles using the BackTracker
/**
 * This is a class to make matching reco objects to MCParticles
 * easier. The cheat::BackTracker class is used to do the matching.
 *
 * fcl configuration options:
 *
 * HitTag: Input tag of the hit collection used to make the reco 
 *   objects.
 *
 */
class mctrue::BackTrackMatcherAlg
{
  public:
    BackTrackMatcherAlg(fhicl::ParameterSet const& pset);
    void reconfigure(fhicl::ParameterSet const& pset);

      const simb::MCParticle* getBestMatch(
                art::Ptr<recob::Track> const track,
                art::Event const& event);

//    /// Find best match MCParticle given track based on N hits
//    /**
//     * Uses backtracker to find the MCParticle that best matches
//     * the given track based on hit purity i.e. the MCParticle
//     * that made the most hits that match to the track.
//     *
//     * allHits are all hits in the event.
//     */
//    const simb::MCParticle* getBestMatch(
//                art::Ptr<recob::Track> const track,
//                std::vector< art::Ptr<recob::Hit> > const& allHits);
//
//    /// Find best match MCParticle given track based on N hits, with diagnositc variables
//    /**
//     * Uses backtracker to find the MCParticle that best matches
//     * the given track based on hit purity i.e. the MCParticle
//     * that made the most hits that match to the track.
//     *
//     * allHits are all hits in the event.
//     *
//     * hitPur: the hit purity of the track-mcparticle match 
//     *   will be put in the referenced variable. It is the 
//     *   number of hits matched to both the track and mcparticle 
//     *   divided by the number of hits matched to the track.
//     *
//     * hitEnergyPur: the charge purity of hits of the 
//     *   track-mcparticle match will be put in the referenced 
//     *   variable. It is the charge of all of the hits matched
//     *   to both the track and mcparticle divided by the charge
//     *   of the hits matched to the track.
//     *
//     * hitEff: the hit efficiency of the track-mcparticle match 
//     *   will be put in the referenced variable. It is the 
//     *   hits matched to both the track and mcparticle divided 
//     *   by all mcparticle matched hits
//     *
//     * hitEnergyEff: the charge efficiency of hits of the 
//     *   track-mcparticle match will be put in the referenced 
//     *   variable. It is the charge of all of the
//     *   hits matched to both the track and mcparticle divided 
//     *   by the charge of all the mcparticle matched hits
//     */
//    art::Ptr<simb::MCParticle> getBestMatch(
//                art::Ptr<recob::Track> const track,
//                std::vector< art::Ptr<recob::Hit> > const& allHits,
//                float& hitPur, float& hitEnergyPur,
//                float& hitEff, float& hitEnergyEff);
//
//    /// Find best match MCParticle given track based on charge of hits
//    /**
//     * Uses backtracker to find the MCParticle that best matches
//     * the given track based on hit charge purity i.e. the MCParticle
//     * that made the most hit charge that match to the track.
//     *
//     * allHits are all hits in the event.
//     */
//    art::Ptr<simb::MCParticle> getBestMatchCharge(
//                art::Ptr<recob::Track> const track,
//                std::vector< art::Ptr<recob::Hit> > const& allHits);
//
//    /// Find best match MCParticle given track based on charge of hits, with diagnositc variables
//    /**
//     * Uses backtracker to find the MCParticle that best matches
//     * the given track based on hit charge purity i.e. the MCParticle
//     * that made the most hit charge that match to the track.
//     *
//     * allHits are all hits in the event.
//     *
//     * hitPur: the hit purity of the track-mcparticle match 
//     *   will be put in the referenced variable. It is the 
//     *   number of hits matched to both the track and mcparticle 
//     *   divided by the number of hits matched to the track.
//     *
//     * hitEnergyPur: the charge purity of hits of the 
//     *   track-mcparticle match will be put in the referenced 
//     *   variable. It is the charge of all of the hits matched
//     *   to both the track and mcparticle divided by the charge
//     *   of the hits matched to the track.
//     *
//     * hitEff: the hit efficiency of the track-mcparticle match 
//     *   will be put in the referenced variable. It is the 
//     *   hits matched to both the track and mcparticle divided 
//     *   by all mcparticle matched hits
//     *
//     * hitEnergyEff: the charge efficiency of hits of the 
//     *   track-mcparticle match will be put in the referenced 
//     *   variable. It is the charge of all of the
//     *   hits matched to both the track and mcparticle divided 
//     *   by the charge of all the mcparticle matched hits
//     */
//    art::Ptr<simb::MCParticle> getBestMatchCharge(
//                art::Ptr<recob::Track> const track,
//                std::vector< art::Ptr<recob::Hit> > const& allHits,
//                float& hitPur, float& hitEnergyPur,
//                float& hitEff, float& hitEnergyEff);
//
//    /// Find tracks that match the given MCParticle
//    /**
//     * Uses backtracker to find tracks that match the MCParticle
//     * with hit purity greater than minHitPur and hit energy 
//     * purity greater than minHitEnergyPur. The hit purity is
//     * the number of hits matched to both the track and mcparticle
//     * divided by the number of hits matched to the track. The
//     * hit energy purity is the charge of the hits matched to both
//     * the track and mcparticle divided by the charge of the hits
//     * matched to the track.
//     *
//     * allHits are all hits in the event.
//     *
//     * hitEff: the hit efficiency of the track-mcparticle match 
//     *   will be put in the referenced variable. It is the 
//     *   hits matched to both the track and mcparticle divided 
//     *   by all mcparticle matched hits
//     *
//     * hitEnergyEff: the charge efficiency of hits of the 
//     *   track-mcparticle match will be put in the referenced 
//     *   variable. It is the charge of all of the
//     *   hits matched to both the track and mcparticle divided 
//     *   by the charge of all the mcparticle matched hits
//     */
//    std::vector<art::Ptr<recob::Track> > getMatchedTracks(
//                simb::MCParticle const& mcparticle, 
//                std::vector< art::Ptr<recob::Track> > const& alltracks, 
//                std::vector< art::Ptr<recob::Hit> > const& allHits,
//                float minHitPur, float minHitEnergyPur,
//                float& hitEff, float& hitEnergyEff);
//
//    /// Find clusters that match the given MCParticle
//    /**
//     * Uses backtracker to find clusters that match the MCParticle
//     * with hit purity greater than minHitPur and hit energy 
//     * purity greater than minHitEnergyPur. The hit purity is
//     * the number of hits matched to both the cluster and mcparticle
//     * divided by the number of hits matched to the cluster. The
//     * hit energy purity is the charge of the hits matched to both
//     * the cluster and mcparticle divided by the charge of the hits
//     * matched to the cluster.
//     *
//     * allHits are all hits in the event.
//     *
//     * hitEff: the hit efficiency of the track-mcparticle match 
//     *   will be put in the referenced variable. It is the 
//     *   hits matched to both the cluster and mcparticle divided 
//     *   by all mcparticle matched hits
//     *
//     * hitEnergyEff: the charge efficiency of hits of the 
//     *   track-mcparticle match will be put in the referenced 
//     *   variable. It is the charge of all of the
//     *   hits matched to both the cluster and mcparticle divided 
//     *   by the charge of all the mcparticle matched hits
//     */
//    std::vector<art::Ptr<recob::Cluster> > getMatchedClusters(
//                simb::MCParticle const& mcparticle, 
//                std::vector< art::Ptr<recob::Cluster> > const& allclusters,
//                std::vector< art::Ptr<recob::Hit> > const& allHits,
//                float minHitPur, float minHitEnergyPur,
//                float& hitEff, float& hitEnergyEff);
//
//    /// Find hits that match the given mcparticle
//    /**
//     * Uses backtracker to find hits that match the MCParticle
//     */
//    std::vector<art::Ptr<recob::Hit> > getMatchedHits(
//                simb::MCParticle const& mcparticle,
//                std::vector< art::Ptr<recob::Hit> > const& allHits);
//
//    /// Find the fraction of the mcparticle's charge reconstructed in hits
//    /**
//     * Uses backtracker to find IDEs and hits that match the 
//     * MCParticle, and find the fraction of the IDE charge that 
//     * gets reconstructed as hit charge.
//     */
//    float getHitChargeEfficiency(
//                simb::MCParticle const& mcparticle, 
//                std::vector< art::Ptr<recob::Hit> > const& allHits);
//
    /// Return (efficiency,purity) of hits for mcparticle
    /**
     * Returns a pair (efficiency,purity) of the given charge 
     * in the given hits matching to the given mcparticle.
     * allHits are all reconstructed hits and used as the 
     * denominator of the efficiency.
     */
    std::pair<double,double> getHitEffPur(
                simb::MCParticle const& mcparticle,
    		    std::vector< art::Ptr<recob::Hit> > const& hits,
    		    std::vector< art::Ptr<recob::Hit> > const& allHits);

    /// Return charge (efficiency,purity) of hits for mcparticle
    /**
     * Returns a pair (efficiency,purity) of the given hits 
     * matching to the given mcparticle. allHits are all 
     * reconstructed hits and used as the denominator of the
     * efficiency.
     */
    std::pair<double,double> getHitChargeEffPur(
                simb::MCParticle const& mcparticle,
    		    std::vector< art::Ptr<recob::Hit> > const& hits,
    		    std::vector< art::Ptr<recob::Hit> > const& allHits);

  private:

    /// Find best match MCParticle given these hits, with diagnositc variables
    /**
     * Uses backtracker to find the MCParticle that best matches
     * theseHits based on hit purity i.e. the MCParticle
     * that made the most of these hits or hit charge purity i.e. 
     * the MCParticle that left the most charge in these hits.
     *
     * allHits are all hits in the event.
     *
     * hitPur: the hit purity of the track-mcparticle match 
     *   will be put in the referenced variable. It is the 
     *   number of theseHits matched to the mcparticle 
     *   divided by the number of theseHits.
     *
     * hitEnergyPur: the charge purity of hits of the 
     *   track-mcparticle match will be put in the referenced 
     *   variable. It is the charge of all of theseHits matched
     *   to the mcparticle divided by the charge of all of 
     *   theseHits.
     *
     * hitEff: the hit efficiency of the track-mcparticle match 
     *   will be put in the referenced variable. It is the 
     *   hits matched to both the track and mcparticle divided 
     *   by all mcparticle matched hits
     *
     * hitEnergyEff: the charge efficiency of hits of the 
     *   track-mcparticle match will be put in the referenced 
     *   variable. It is the charge of all of the
     *   hits matched to both the track and mcparticle divided 
     *   by the charge of all the mcparticle matched hits
     *
     * charge: if true, match based on charge, if false, match
     *   based on number of hits.
     */
    const simb::MCParticle* getBestMatchTheseHits(
            std::vector< art::Ptr<recob::Hit> > const& theseHits,
            bool charge);

    art::ServiceHandle<cheat::BackTracker> fBT; ///< the back tracker service
    art::InputTag fHitTag;
    
};

#endif
