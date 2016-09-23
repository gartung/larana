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
#include <set>

truthmatching::BackTrackMatcherAlg::BackTrackMatcherAlg(fhicl::ParameterSet const& pset): fAllHits(), fIsSetup(false), fBT()
{
  return;
} //end constructor

void 
truthmatching::BackTrackMatcherAlg::reconfigure(fhicl::ParameterSet const& pset)
{
  return;
}

void
truthmatching::BackTrackMatcherAlg::setup(std::vector< art::Ptr<recob::Hit> > const& allhits)
{
  fIsSetup = true;
  fAllHits = allhits;
}

art::Ptr<recob::Track>
truthmatching::BackTrackMatcherAlg::getBestMatch(simb::MCParticle const& mcparticle, 
                    std::vector< art::Ptr<recob::Track> > const& alltracks)
{
  return art::Ptr<recob::Track>();
}

art::Ptr<recob::Track>
truthmatching::BackTrackMatcherAlg::getBestMatch(simb::MCParticle const& mcparticle, 
                    std::vector< art::Ptr<recob::Track> > const& alltracks, 
                    float& hitEff, float& hitPur, // result diagnostic values
                    float& hitEnergyEff, float& hitEnergyPur) // result diagnostic values
{
  return art::Ptr<recob::Track>();
}


std::vector<art::Ptr<recob::Cluster> >
truthmatching::BackTrackMatcherAlg::getMatchedClusters(simb::MCParticle const& mcparticle, 
                std::vector< art::Ptr<recob::Cluster> > const& allclusters,
                std::vector<float>& hitEff, std::vector<float>& hitPur, // result diagnostic values
                std::vector<float>& hitEnergyEff, std::vector<float>& hitEnergyPur) // result diagnostic values
{
  return std::vector<art::Ptr<recob::Cluster> >();
}

std::vector<art::Ptr<recob::Hit> >
truthmatching::BackTrackMatcherAlg::getMatchedHits(simb::MCParticle const& mcparticle)
{
  return std::vector<art::Ptr<recob::Hit> >();
}

float 
truthmatching::BackTrackMatcherAlg::getHitEnergyEfficiency(simb::MCParticle const& mcparticle)
{
  return -1.;
}

std::pair<double,double>
truthmatching::BackTrackMatcherAlg::getHitEffPur(simb::MCParticle const& mcparticle,
				std::vector< art::Ptr<recob::Hit> > const& hits)
{
    const geo::View_t view = hits[0]->View();
    const int trackID = mcparticle.TrackId();
    std::set<int> id;
    id.insert(trackID);
    
    // use the cheat::BackTracker to find purity and efficiency for these hits
    double purity     = fBT->HitCollectionPurity(id, hits);
    double efficiency = fBT->HitCollectionEfficiency(id, hits, fAllHits, view);
    return std::make_pair(purity,efficiency);
}

std::pair<double,double>
truthmatching::BackTrackMatcherAlg::getHitChargeEffPur(
                simb::MCParticle const& mcparticle,
				std::vector< art::Ptr<recob::Hit> > const& hits)
{
    const geo::View_t view = hits[0]->View();
    const int trackID = mcparticle.TrackId();
    std::set<int> id;
    id.insert(trackID);

    // use the cheat::BackTracker to find purity and efficiency for these hits
    double purity     = fBT->HitChargeCollectionPurity(id, hits);
    double efficiency = fBT->HitChargeCollectionEfficiency(id, hits, fAllHits, view);
    return std::make_pair(purity,efficiency);
}
