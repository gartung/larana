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

    /// Find best match MCParticle given recoObj based on N hits
    /**
     * Uses backtracker to find the MCParticle that best matches
     * the given reco object based on hit purity i.e. the MCParticle
     * that made the most hits that match to the recoobj.
     */
      template<typename T> inline
      const simb::MCParticle* getBestMatch(
                art::Ptr<T> const recoObj,
                art::Event const& event);

    /// Find best match MCParticle given recoObj based on N hits, with diagnostic values
    /**
     * Uses backtracker to find the MCParticle that best matches
     * the given reco object based on hit purity i.e. the MCParticle
     * that made the most hits that match to the recoObj.
     *
     * hitPur: the hit purity of the recoObj-mcparticle match 
     *   will be put in the referenced variable. It is the 
     *   number of hits matched to both the recoOjb and mcparticle 
     *   divided by the number of hits matched to the recoOjb.
     *
     * hitEnergyPur: the charge purity of hits of the 
     *   recoObj-mcparticle match will be put in the referenced 
     *   variable. It is the charge of all of the hits matched
     *   to both the recoObj and mcparticle divided by the charge
     *   of the hits matched to the recoObj.
     *
     * hitEff: the hit efficiency of the recoObj-mcparticle match 
     *   will be put in the referenced variable. It is the 
     *   hits matched to both the recoObj and mcparticle divided 
     *   by all mcparticle matched hits
     *
     * hitEnergyEff: the charge efficiency of hits of the 
     *   recoObj-mcparticle match will be put in the referenced 
     *   variable. It is the charge of all of the
     *   hits matched to both the recoObj and mcparticle divided 
     *   by the charge of all the mcparticle matched hits
     */
      template<typename T> inline
      const simb::MCParticle* getBestMatch(
                art::Ptr<T> const recoObj,
                art::Event const& event,
                float& hitPur, float& hitEnergyPur,
                float& hitEff, float& hitEnergyEff);

    /// Find best match MCParticle given recoObj based on charge of hits
    /**
     * Uses backtracker to find the MCParticle that best matches
     * the given reco object based on charge purity i.e. the MCParticle
     * that made the most hit charge that match to the recoobj.
     *
     */
      template<typename T> inline
      const simb::MCParticle* getBestMatchCharge(
                art::Ptr<T> const recoObj,
                art::Event const& event);

    /// Find best match MCParticle given recoObj based on charge of hits, with diagnostic values
    /**
     * Uses backtracker to find the MCParticle that best matches
     * the given reco object based on charge purity i.e. the MCParticle
     * that made the most hit charge that match to the recoobj.
     *
     * hitPur: the hit purity of the recoObj-mcparticle match 
     *   will be put in the referenced variable. It is the 
     *   number of hits matched to both the recoOjb and mcparticle 
     *   divided by the number of hits matched to the recoOjb.
     *
     * hitEnergyPur: the charge purity of hits of the 
     *   recoObj-mcparticle match will be put in the referenced 
     *   variable. It is the charge of all of the hits matched
     *   to both the recoObj and mcparticle divided by the charge
     *   of the hits matched to the recoObj.
     *
     * hitEff: the hit efficiency of the recoObj-mcparticle match 
     *   will be put in the referenced variable. It is the 
     *   hits matched to both the recoObj and mcparticle divided 
     *   by all mcparticle matched hits
     *
     * hitEnergyEff: the charge efficiency of hits of the 
     *   recoObj-mcparticle match will be put in the referenced 
     *   variable. It is the charge of all of the
     *   hits matched to both the recoObj and mcparticle divided 
     *   by the charge of all the mcparticle matched hits
     */
      template<typename T> inline
      const simb::MCParticle* getBestMatchCharge(
                art::Ptr<T> const recoObj,
                art::Event const& event,
                float& hitPur, float& hitEnergyPur,
                float& hitEff, float& hitEnergyEff);

    /// Find recoObjs that match the given MCParticle
    /**
     * Uses backtracker to find tracks that match the MCParticle
     * with hit purity greater than minHitPur and hit energy 
     * purity greater than minHitEnergyPur. The hit purity is
     * the number of hits matched to both the track and mcparticle
     * divided by the number of hits matched to the track. The
     * hit energy purity is the charge of the hits matched to both
     * the track and mcparticle divided by the charge of the hits
     * matched to the track.
     *
     * allHits are all hits in the event.
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
     */
      template<typename T> inline
      std::vector<art::Ptr<T> > getMatched(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> > const& recoObjs, 
                art::Event const& event,
                float minHitPur, float minHitEnergyPur);

      template<typename T> inline
      std::vector<art::Ptr<T> > getMatched(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> > const& recoObjs, 
                art::Event const& event,
                float minHitPur, float minHitEnergyPur,
                float& hitEff, float& hitEnergyEff);

      /// Sorts reco objects by distance to MCParticle start point
      /**
        * If you feed a vector of art::Ptr's of reco objects in 
        * here, it will sort them by nearest hit to the MCParticle
        * start point, nearest object first. This uses BackTracker
        * to find the object's hits and then finds which simIDEs
        * those hits came from. The nearest simIDE to the MCParticle
        * start point is then found. The objects are then sorted
        * by this distance, closest to farthest. The sorted vector
        * is returned.
        *
        * This is a fairly slow algorithm because it has to compute
        * the distance to every hit.
        */
      template<typename T> inline
      std::vector<art::Ptr<T> > sortByDistToMCParticleStart(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> >& recoObjs,
                art::Event const& event);

      /// Sorts reco objects by distance to MCParticle end point
      /**
        * If you feed a vector of art::Ptr's of reco objects in 
        * here, it will sort them by nearest hit to the MCParticle
        * end point, nearest object first. This uses BackTracker
        * to find the object's hits and then finds which simIDEs
        * those hits came from. The nearest simIDE to the MCParticle
        * end point is then found. The objects are then sorted
        * by this distance, closest to farthest. The sorted vector
        * is returned.
        *
        * This is a fairly slow algorithm because it has to compute
        * the distance to every hit.
        */
      template<typename T> inline
      std::vector<art::Ptr<T> > sortByDistToMCParticleEnd(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> >& recoObjs,
                art::Event const& event);

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

    /// Find indices of reco objs that match mcparticle, with diagnositc variables
    /**
     * Uses backtracker to find tracks that match  MCParticle that best matches
     * with minimum hit and charge purity.
     */
    const std::vector<size_t> getMatchedSetsOfHits(
                simb::MCParticle const& mcparticle, 
                art::FindManyP<recob::Hit>& fmh,
                float minHitPur, float minHitEnergyPur,
                float& hitEff, float& hitEnergyEff);

    /// Find sorted indices of sets of hits, based on distance to point
    /**
     * Returns a vector of indices, sorting the input sets of hits by
     * distance to the given point. Uses BackTracker to find the
     * true positions of the hits.
     */
    const std::vector<size_t> sortObjsByDistance(
                const TVector3 & point, 
                const art::FindManyP<recob::Hit> & fmh);

    //////////////////////
    // Member Variables //
    //////////////////////

    art::ServiceHandle<cheat::BackTracker> fBT; ///< the back tracker service
    art::InputTag fHitTag; ///< input tag of hits used by this algo
    
};

//////////////////////////
// Template Definitions //
//////////////////////////

template<typename T> inline
const simb::MCParticle* 
mctrue::BackTrackMatcherAlg::getBestMatch(
          art::Ptr<T> const recoObj,
          art::Event const& event)
{
  float hitPur=0; 
  float hitEnergyPur=0;
  float hitEff=0; 
  float hitEnergyEff=0;
  return getBestMatch(recoObj,event,hitPur,hitEnergyPur,hitEff,hitEnergyEff);
}

template<typename T> inline
const simb::MCParticle* 
mctrue::BackTrackMatcherAlg::getBestMatch(
          art::Ptr<T> const recoObj,
          art::Event const& event,
          float& hitPur, float& hitEnergyPur,
          float& hitEff, float& hitEnergyEff)
{
  art::FindManyP<recob::Hit> fmh({recoObj}, event, fHitTag);
  const std::vector<art::Ptr<recob::Hit>> & hits = fmh.at(0);
  return getBestMatchTheseHits(hits,false);
}

template<typename T> inline
const simb::MCParticle* 
mctrue::BackTrackMatcherAlg::getBestMatchCharge(
          art::Ptr<T> const recoObj,
          art::Event const& event)
{
  float hitPur=0; 
  float hitEnergyPur=0;
  float hitEff=0; 
  float hitEnergyEff=0;
  return getBestMatchCharge(recoObj,event,hitPur,hitEnergyPur,hitEff,hitEnergyEff);
}

template<typename T> inline
const simb::MCParticle* 
mctrue::BackTrackMatcherAlg::getBestMatchCharge(
          art::Ptr<T> const recoObj,
          art::Event const& event,
          float& hitPur, float& hitEnergyPur,
          float& hitEff, float& hitEnergyEff)
{
  art::FindManyP<recob::Hit> fmh({recoObj}, event, fHitTag);
  const std::vector<art::Ptr<recob::Hit>> & hits = fmh.at(0);
  return getBestMatchTheseHits(hits,true);
}

template<typename T> inline
std::vector<art::Ptr<T> > 
mctrue::BackTrackMatcherAlg::getMatched(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> > const& recoObjs, 
                art::Event const& event,
                float minHitPur, float minHitEnergyPur)
{
  float hitEff=0; 
  float hitEnergyEff=0;
  return getMatched(mcparticle,recoObjs,event,minHitPur,minHitEnergyPur,hitEff,hitEnergyEff);
}

template<typename T> inline
std::vector<art::Ptr<T> > 
mctrue::BackTrackMatcherAlg::getMatched(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> > const& recoObjs, 
                art::Event const& event,
                float minHitPur, float minHitEnergyPur,
                float& hitEff, float& hitEnergyEff)
{
  std::vector<art::Ptr<recob::Track>> result;
  art::FindManyP<recob::Hit> fmh(recoObjs, event, fHitTag);
  const std::vector<size_t> matchedIs = getMatchedSetsOfHits(mcparticle,fmh,minHitPur,
                        minHitEnergyPur,hitEff,hitEnergyEff);
  for (auto matchedI: matchedIs)
  {
    result.push_back(recoObjs.at(matchedI));
  }
  return result;
}

//// Find recoObjs that match the given MCParticle (specialization for hits)
//template<> inline
//std::vector<art::Ptr<recob::Hit> > 
//mctrue::BackTrackMatcherAlg::getMatched(
//          simb::MCParticle const& mcparticle, 
//          std::vector< art::Ptr<recob::Hit> > const& hits, 
//          art::Event const& event,
//          float minHitPur, float minHitEnergyPur,
//          float& hitEff, float& hitEnergyEff)
//{
//  hitEff=-1; 
//  hitEnergyEff=-1;
//  const auto matchedHitsVec = bt->TrackIDsToHits(hits,{mcparticle.TrackId()});
//  return matchedHitsVec.at(0);
//}

template<typename T> inline
std::vector<art::Ptr<T> >
mctrue::BackTrackMatcherAlg::sortByDistToMCParticleStart(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> >& recoObjs,
                art::Event const& event)
{
  art::FindManyP<recob::Hit> fmh(recoObjs, event, fHitTag);
  const std::vector<size_t> indices = sortObjsByDistance(mcparticle.Position().Vect(),recoObjs,fmh);
  std::vector<art::Ptr<T>> result;
  for (const auto& index: indices)
  {
    result.push_back(recoObjs[index]);
  }
  return result;
}

template<typename T> inline
std::vector<art::Ptr<T> >
mctrue::BackTrackMatcherAlg::sortByDistToMCParticleEnd(
                simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<T> >& recoObjs,
                art::Event const& event)
{
  art::FindManyP<recob::Hit> fmh(recoObjs, event, fHitTag);
  const std::vector<size_t> indices = sortObjsByDistance(mcparticle.EndPosition().Vect(),fmh);
  std::vector<art::Ptr<T>> result;
  for (const auto& index: indices)
  {
    result.push_back(recoObjs[index]);
  }
  return result;
}

#endif
