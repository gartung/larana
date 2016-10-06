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
#include <algorithm>
#include <iostream>

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
    return std::make_pair(efficiency,purity);
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
    return std::make_pair(efficiency,purity);
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
  std::cout << "BackTrackMatcherAlg::getMatchedSetsOfHits called\n";
  hitEff = -1;
  hitEnergyEff = -1;
  std::vector<size_t> result;
  const std::vector< art::Ptr<recob::Hit> > allHits;
  for (unsigned iTrack=0; iTrack<fmh.size(); iTrack++)
  {
    std::cout << "iTrack: "<< iTrack << "\n";
    const std::vector< art::Ptr<recob::Hit> > thisTrackHits = fmh.at(iTrack);
    std::cout << "nhits: "<< thisTrackHits.size() << "\n";
    auto hiteffpur = getHitEffPur(mcparticle,thisTrackHits,allHits);
    std::cout << "hiteff: "<< hiteffpur.first << " hitpur: "<< hiteffpur.second << "\n";
    if(hiteffpur.second < minHitPur)
    {
      continue;
    }
    auto chargeeffpur = getHitChargeEffPur(mcparticle,thisTrackHits,allHits);
    std::cout << "chargeeff: "<< chargeeffpur.first << " chargepur: "<< chargeeffpur.second << "\n";
    if(chargeeffpur.second < minHitEnergyPur)
    {
      continue;
    }
    result.push_back(iTrack);
  }
  return result;
}

const std::vector<size_t> 
mctrue::BackTrackMatcherAlg::sortObjsByDistance(
            const TVector3 & point, 
            const art::FindManyP<recob::Hit> & fmh)
{
  std::vector<size_t> result;
  std::vector<double> objDistances;
  for (size_t iObj=0; iObj<fmh.size(); iObj++)
  {
    const std::vector< art::Ptr<recob::Hit> > objHits = fmh.at(iObj);
    double closestDistance = 1.0e9;
    for (const auto & hit: objHits)
    {
      std::vector<double> hitPointVect = fBT->HitToXYZ(hit);
      if(hitPointVect.size() != 3)  
      {
        throw cet::exception("IvalidSize","HitToXYZ returned a vector with size != 3, size: ");
      }
      const TVector3 hitPoint(hitPointVect[0],hitPointVect[1],hitPointVect[2]);
      const double hitDistance = (hitPoint - point).Mag();
      if (hitDistance < closestDistance)
      {
        closestDistance = hitDistance;
      }
    } // for hit
    objDistances.push_back(closestDistance);
    result.push_back(iObj);
  } // for iObj
  std::sort(result.begin(),result.end(),
            [&objDistances](size_t a, size_t b){
                return (objDistances[a]<objDistances[b]);
            }
  );
  return result;
}
