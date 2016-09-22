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

truthmatching::BackTrackMatcherAlg::BackTrackMatcherAlg(fhicl::ParameterSet const& pset): fAllHits(), fIsSetup(false)
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
truthmatching::BackTrackMatcherAlg::getBestMatch(simb::MCParticle & mcparticle, 
                    std::vector< art::Ptr<recob::Track> > const& alltracks)
{
  return art::Ptr<recob::Track>();
}

art::Ptr<recob::Track>
truthmatching::BackTrackMatcherAlg::getBestMatch(simb::MCParticle & mcparticle, 
                    std::vector< art::Ptr<recob::Track> > const& alltracks, 
                    float& hitEff, float& hitPur, // result diagnostic values
                    float& hitEnergyEff, float& hitEnergyPur) // result diagnostic values
{
  return art::Ptr<recob::Track>();
}


std::vector<art::Ptr<recob::Cluster> >
truthmatching::BackTrackMatcherAlg::getMatchedClusters(simb::MCParticle & mcparticle, 
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
truthmatching::BackTrackMatcherAlg::getHitEnergyEfficiency(simb::MCParticle & mcparticle)
{
  return -1.;
}

