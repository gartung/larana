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

namespace truthmatching {
  class BackTrackMatcherAlg;
}

class truthmatching::BackTrackMatcherAlg
{
  public:
    BackTrackMatcherAlg(fhicl::ParameterSet const& pset);
    void reconfigure(fhicl::ParameterSet const& pset);

    void setup(std::vector< art::Ptr<recob::Hit> > const& allhits); // setup every event before calling other methods

    art::Ptr<recob::Track> getBestMatch(simb::MCParticle const& mcparticle, 
                  std::vector< art::Ptr<recob::Track> > const& alltracks);

    art::Ptr<recob::Track> getBestMatch(simb::MCParticle const& mcparticle, 
                        std::vector< art::Ptr<recob::Track> > const& alltracks, 
                        float& hitEff, float& hitPur, // result diagnostic values
                        float& hitEnergyEff, float& hitEnergyPur); // result diagnostic values

    std::vector<art::Ptr<recob::Cluster> > getMatchedClusters(simb::MCParticle const& mcparticle, 
              std::vector< art::Ptr<recob::Cluster> > const& allclusters,
              std::vector<float>& hitEff, std::vector<float>& hitPur, // result diagnostic values
              std::vector<float>& hitEnergyEff, std::vector<float>& hitEnergyPur); // result diagnostic values

    std::vector<art::Ptr<recob::Hit> > getMatchedHits(simb::MCParticle const& mcparticle);

    float getHitEnergyEfficiency(simb::MCParticle const& mcparticle);

  private:

    std::pair<double,double> getHitEffPur(simb::MCParticle const& mcparticle,
    				std::vector< art::Ptr<recob::Hit> > const& hits);

    std::pair<double,double> getHitChargeEffPur(
                simb::MCParticle const& mcparticle,
				std::vector< art::Ptr<recob::Hit> > const& hits);
    


    std::vector< art::Ptr<recob::Hit> > fAllHits;
    bool fIsSetup;
    art::ServiceHandle<cheat::BackTracker> fBT; ///< the back tracker service
    
};

#endif
