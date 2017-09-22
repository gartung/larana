////////////////////////////////////////////////////////////////////////
// Class:       CosmicTaggerEfficiency
// Plugin Type: analyzer (art v2_07_03)
// File:        CosmicTaggerEfficiency_module.cc
//
// Generated at Fri Sep 22 09:50:29 2017 by Tingjun Yang using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larsim/MCCheater/BackTracker.h"

#include "TEfficiency.h"
#include "TH1F.h"

namespace cosmic {
  class CosmicTaggerEfficiency;
}


class cosmic::CosmicTaggerEfficiency : public art::EDAnalyzer {
public:
  explicit CosmicTaggerEfficiency(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicTaggerEfficiency(CosmicTaggerEfficiency const &) = delete;
  CosmicTaggerEfficiency(CosmicTaggerEfficiency &&) = delete;
  CosmicTaggerEfficiency & operator = (CosmicTaggerEfficiency const &) = delete;
  CosmicTaggerEfficiency & operator = (CosmicTaggerEfficiency &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  int GetOrigin(std::vector<art::Ptr<recob::Hit>> allHits);

private:

  // Declare member data here.
  std::vector<art::InputTag> fCosmicProducerLabels;     ///< List of cosmic tagger producers
  art::InputTag              fTrackProducerLabel;       ///< Track producer
  std::vector<float>         fCosmicTagThresholds;      ///< Thresholds for tagging

  TEfficiency* eff_cosmic = 0;
  TEfficiency* eff_neutrino = 0;

  TH1F *h_cosmic;
  TH1F *h_cosmic_tag;
  TH1F *h_neutrino;
  TH1F *h_neutrino_tag;

};


cosmic::CosmicTaggerEfficiency::CosmicTaggerEfficiency(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCosmicProducerLabels(p.get<std::vector<art::InputTag>>("CosmicProducerLabels")),
  fTrackProducerLabel  (p.get<art::InputTag>             ("TrackProducerLabel")),
  fCosmicTagThresholds (p.get<std::vector<float>>        ("CosmicTagThresholds"))
{
}

void cosmic::CosmicTaggerEfficiency::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  std::map<unsigned, int> originmap;
  std::map<unsigned, bool> tagmap;

  // Tracks
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel(fTrackProducerLabel,trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(trackHandle, evt, fTrackProducerLabel);

  for (size_t i = 0; i<tracks.size(); ++i){
    art::Ptr<recob::Track> track = tracks[i];
    std::vector< art::Ptr<recob::Hit> > allHits = fmh.at(i);
    int origin = GetOrigin(allHits);
    bool tagged = false;
    if (origin == 1) tagged = true; //neutrino
    for (size_t j = 0; j<fCosmicProducerLabels.size(); ++j){
      art::FindManyP<anab::CosmicTag> fmc(trackHandle, evt, fCosmicProducerLabels[j]);
      std::vector<art::Ptr<anab::CosmicTag>> cosmictags = fmc.at(i);
      if (cosmictags.size()){
        if (cosmictags[0]->CosmicScore()>fCosmicTagThresholds[j]){
          if (origin == 2) tagged = true;  //cosmics
          if (origin == 1) tagged = false; //neutrinos
        }
      }
    }
    if (origin == 1){
      h_neutrino->Fill(track->Length());
      if (tagged) h_neutrino_tag->Fill(track->Length());
    }
    if (origin == 2){
      h_cosmic->Fill(track->Length());
      if (tagged) h_cosmic_tag->Fill(track->Length());
    }
  }

}

void cosmic::CosmicTaggerEfficiency::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  h_cosmic = tfs->make<TH1F>("h_cosmic","Cosmic Tracks;Track length (cm);Tracks", 200,0,1000);
  h_cosmic_tag = tfs->make<TH1F>("h_cosmic_tag","Tagged Cosmic Tracks;Track length (cm);Tracks", 200,0,1000);
  h_neutrino = tfs->make<TH1F>("h_neutrino","Neutrino Tracks;Track length (cm);Tracks", 200,0,1000);
  h_neutrino_tag = tfs->make<TH1F>("h_neutrino_tag","Tagged Neutrino Tracks;Track length (cm);Tracks", 200,0,1000);
}

void cosmic::CosmicTaggerEfficiency::endJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  if(TEfficiency::CheckConsistency(*h_cosmic_tag,*h_cosmic)){
    eff_cosmic = tfs->make<TEfficiency>(*h_cosmic_tag,*h_cosmic);
    eff_cosmic->Write("eff_cosmic");
  }
  if(TEfficiency::CheckConsistency(*h_neutrino_tag,*h_neutrino)){
    eff_neutrino = tfs->make<TEfficiency>(*h_neutrino_tag,*h_neutrino);
    eff_neutrino->Write("eff_neutrino");
  }

}

int cosmic::CosmicTaggerEfficiency::GetOrigin(std::vector<art::Ptr<recob::Hit>> allHits){
  
  art::ServiceHandle<cheat::BackTracker> bt;

  std::map<int,double> trkide;
  for(size_t h = 0; h < allHits.size(); ++h){
    art::Ptr<recob::Hit> hit = allHits[h];
    //std::vector<sim::IDE> ides;
    std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
    
    for(size_t e = 0; e < TrackIDs.size(); ++e){
      trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
    }
  }
  // Work out which IDE despoited the most charge in the hit if there was more than one.
  double maxe = -1;
  double tote = 0;
  int Trackid = INT_MIN;
  
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      //if(pfPartIdx < max_pfparticles) origin=ii->first;
      Trackid=ii->first;
    }
  }
  
  if (Trackid!=INT_MIN){
    const simb::MCParticle* particle=bt->TrackIDToParticle(Trackid);	    
    if(particle){
      return bt->TrackIDToMCTruth(Trackid)->Origin();
    }
    else return 0;
  }
  else return 0;
}

DEFINE_ART_MODULE(cosmic::CosmicTaggerEfficiency)
