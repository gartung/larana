/*!
 * Title:   BackTrack Matcher Algorithim Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for matching MCParticles to reconstructed hits, clusters, and tracks
 * Input:       simb::MCParticle and (vector<recob::Hit> or vector<recob::Cluster> or vector<recob::Track>)
 * Output:      vector<recob::Hit> or vector<recob::Cluster> or recob::Track and matching info
*/

#include "BackTrackMatcherAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include <set>

mctrue::BackTrackMatcherAlg::BackTrackMatcherAlg(fhicl::ParameterSet const& pset): fBT(), fHitTag()
{
  reconfigure(pset);
  return;
} //end constructor

void 
mctrue::BackTrackMatcherAlg::reconfigure(fhicl::ParameterSet const& pset)
{
  fHitTag = pset.get<art::InputTag>("HitTag");
  return;
}

std::pair<double,double>
mctrue::BackTrackMatcherAlg::getHitEffPur(
                    simb::MCParticle const& mcparticle,
    				std::vector< art::Ptr<recob::Hit> > const& hits,
    				std::vector< art::Ptr<recob::Hit> > const& allHits)
{
    const geo::View_t view = hits[0]->View();
    const int trackID = mcparticle.TrackId();
    std::set<int> id;
    id.insert(trackID);
    
    // use the cheat::BackTracker to find purity and efficiency for these hits
    double purity     = fBT->HitCollectionPurity(id, hits);
    double efficiency = fBT->HitCollectionEfficiency(id, hits, allHits, view);
    return std::make_pair(purity,efficiency);
}

std::pair<double,double>
mctrue::BackTrackMatcherAlg::getHitChargeEffPur(
                    simb::MCParticle const& mcparticle,
    				std::vector< art::Ptr<recob::Hit> > const& hits,
    				std::vector< art::Ptr<recob::Hit> > const& allHits)
{
    const geo::View_t view = hits[0]->View();
    const int trackID = mcparticle.TrackId();
    std::set<int> id;
    id.insert(trackID);

    // use the cheat::BackTracker to find purity and efficiency for these hits
    double purity     = fBT->HitChargeCollectionPurity(id, hits);
    double efficiency = fBT->HitChargeCollectionEfficiency(id, hits, allHits, view);
    return std::make_pair(purity,efficiency);
}

const simb::MCParticle*
mctrue::BackTrackMatcherAlg::getBestMatchTheseHits(
            std::vector< art::Ptr<recob::Hit> > const& theseHits,
            bool charge)
{
  std::set<int> trackIDs = fBT->GetSetOfTrackIDs(theseHits);
  int bestTrackID = -1;
  double bestPurity = -1.;
  for(auto trackID: trackIDs)
  {
    std::set<int> id;
    id.insert(trackID);

    double purity = -1.;
    if (charge)
    {
      purity     = fBT->HitCollectionPurity(id, theseHits);
    }
    else
    {
      purity     = fBT->HitCollectionPurity(id, theseHits);
    }

    if(purity > bestPurity)
    {
      bestTrackID = trackID;
      bestPurity = purity;
    }
  } // end loop over trackIDs

  if (bestPurity < 0.) // no MCParticles
  {
    return NULL;
  }
  return fBT->TrackIDToParticle(bestTrackID);
}

const std::vector<size_t> 
mctrue::BackTrackMatcherAlg::getMatchedSetsOfHits(
                simb::MCParticle const& mcparticle, 
                art::FindManyP<recob::Hit>& fmh,
                float minHitPur, float minHitEnergyPur,
                float& hitEff, float& hitEnergyEff)
{
  hitEff = -1;
  hitEnergyEff = -1;
  std::vector<size_t> result;
  const std::vector< art::Ptr<recob::Hit> > allHits;
  for (unsigned iTrack=0; iTrack<fmh.size(); iTrack++)
  {
    const std::vector< art::Ptr<recob::Hit> > thisTrackHits = fmh.at(iTrack);
    auto hiteffpur = getHitEffPur(mcparticle,thisTrackHits,allHits);
    if(hiteffpur.second < minHitPur)
    {
      continue;
    }
    auto chargeeffpur = getHitChargeEffPur(mcparticle,thisTrackHits,allHits);
    if(chargeeffpur.second < minHitEnergyPur)
    {
      continue;
    }
    result.push_back(iTrack);
  }
  return result;
}
