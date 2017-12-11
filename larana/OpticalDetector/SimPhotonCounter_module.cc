// \file SimPhotonCounter.h
// \author Ben Jones, MIT 2010
//
// Module to determine how many phots have been detected at each OpDet
//
// This analyzer takes the SimPhotonsCollection generated by LArG4's sensitive detectors
// and fills up to four trees in the histograms file.  The four trees are:
//
// OpDetEvents       - count how many phots hit the OpDet face / were detected across all OpDet's per event
// OpDets            - count how many phots hit the OpDet face / were detected in each OpDet individually for each event
// AllPhotons      - wavelength information for each phot hitting the OpDet face
// DetectedPhotons - wavelength information for each phot detected
//
// The user may supply a quantum efficiency and sensitive wavelength range for the OpDet's.
// with a QE < 1 and a finite wavelength range, a "detected" phot is one which is
// in the relevant wavelength range and passes the random sampling condition imposed by
// the quantum efficiency of the OpDet
//
// PARAMETERS REQUIRED:
// int32   Verbosity          - whether to write to screen a well as to file. levels 0 to 3 specify different levels of detail to display
// string  InputModule        - the module which produced the SimPhotonsCollection
// bool    MakeAllPhotonsTree - whether to build and store each tree (performance can be enhanced by switching off those not required)
// bool    MakeDetectedPhotonsTree
// bool    MakeOpDetsTree
// bool    MakeOpDetEventsTree
// double  QantumEfficiency   - Quantum efficiency of OpDet
// double  WavelengthCutLow   - Sensitive wavelength range of OpDet
// double  WavelengthCutHigh
#ifndef SimPhotonCounterModule_h
#define SimPhotonCounterModule_h 1

// ROOT includes.
#include "TTree.h"
#include "TFile.h"
#include <Rtypes.h>

// FMWK includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/OpDetResponseInterface.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include <TH1D.h>
#include <TF1.h>
#include <TTree.h>

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace opdet {
  
  class SimPhotonCounter : public art::EDAnalyzer{
    public:
      
      SimPhotonCounter(const fhicl::ParameterSet&);
      virtual ~SimPhotonCounter();
      
      void analyze(art::Event const&);
      
      void beginJob();
      void endJob();
      
    private:
      
      // Trees to output
    
      TTree * fThePhotonTreeAll;
      TTree * fThePhotonTreeDetected;
      TTree * fTheOpDetTree;
      TTree * fTheEventTree;
      
      
      // Parameters to read in
      
      std::string fInputModule;      // Input tag for OpDet collection
      
      int fVerbosity;                // Level of output to write to std::out
      
      bool fMakeDetectedPhotonsTree; //
      bool fMakeAllPhotonsTree;      //
      bool fMakeOpDetsTree;         // Switches to turn on or off each output
      bool fMakeOpDetEventsTree;          //
      
      float fQE;                     // Quantum efficiency of tube
      
      float fWavelengthCutLow;       // Sensitive wavelength range 
      float fWavelengthCutHigh;      // 
      
      TVector3 initialPhotonPosition;
      TVector3 finalPhotonPosition;
      
      
      // Data to store in trees
       
      Float_t fWavelength;
      Float_t fTime;
      Int_t fCount;
      Int_t fCountOpDetAll;
      Int_t fCountOpDetDetected;
      Int_t fCountOpDetReflDetected;
      Float_t fT0_vis;
      
      Int_t fCountEventAll;
      Int_t fCountEventDetected;
      Int_t fCountEventDetectedwithRefl;    
      
      Int_t fEventID;
      Int_t fOpChannel;
      
      //for the analysis tree of the light (gamez)                                                                      
      bool fMakeLightAnalysisTree;
      std::vector<std::vector<std::vector<double> > > fSignals_vuv;
      std::vector<std::vector<std::vector<double> > > fSignals_vis;
      
      TTree * fLightAnalysisTree;
      int fRun, fTrackID, fpdg, fmotherTrackID;
      double fEnergy, fdEdx;
      std::vector<double> fPosition0;
      std::vector<std::vector<double> > fstepPrePositions;
      std::vector<std::vector<double> > fstepPostPositions;
      std::vector<double> fstepPreTimes;
      std::vector<double> fstepPostTimes;
      std::vector<std::vector<double> > fSignalsvuv;
      std::vector<std::vector<double> > fSignalsvis;
      std::string fProcess;
      
      cheat::ParticleInventoryService* pi_serv = nullptr;
  };
}

#endif


namespace opdet {
  
  
  SimPhotonCounter::SimPhotonCounter(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    fVerbosity=                pset.get<int>("Verbosity");
    fInputModule=              pset.get<std::string>("InputModule");
    fMakeAllPhotonsTree=       pset.get<bool>("MakeAllPhotonsTree");
    fMakeDetectedPhotonsTree=  pset.get<bool>("MakeDetectedPhotonsTree");
    fMakeOpDetsTree=           pset.get<bool>("MakeOpDetsTree");
    fMakeOpDetEventsTree=      pset.get<bool>("MakeOpDetEventsTree");
    fMakeLightAnalysisTree=    pset.get<bool>("MakeLightAnalysisTree", false);
    //fQE=                       pset.get<double>("QuantumEfficiency");
    //fWavelengthCutLow=         pset.get<double>("WavelengthCutLow");
    //fWavelengthCutHigh=        pset.get<double>("WavelengthCutHigh");  
  }
  
  
  void SimPhotonCounter::beginJob()
  {
    // Get file service to store trees
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;
    art::ServiceHandle<geo::Geometry> geo;
    
    try { 
      pi_serv = &*(art::ServiceHandle<cheat::ParticleInventoryService>());
    }
    catch (art::Exception const& e) {
      if (e.categoryCode() != art::errors::ServiceNotFound) throw;
      mf::LogError("SimPhotonCounter")
        << "ParticleInventoryService service is not configured!"
        " Please add it in the job configuration."
        " In the meanwhile, some checks to particles will be skipped."
        ;
    }
    
    // Create and assign branch addresses to required tree
    if(fMakeAllPhotonsTree)
    {
      fThePhotonTreeAll = tfs->make<TTree>("AllPhotons","AllPhotons");
      fThePhotonTreeAll->Branch("EventID",     &fEventID,          "EventID/I");
      fThePhotonTreeAll->Branch("Wavelength",  &fWavelength,       "Wavelength/F");
      fThePhotonTreeAll->Branch("OpChannel",       &fOpChannel,            "OpChannel/I");
      fThePhotonTreeAll->Branch("Time",        &fTime,             "Time/F");
    }
    
    if(fMakeDetectedPhotonsTree)
    {
      fThePhotonTreeDetected = tfs->make<TTree>("DetectedPhotons","DetectedPhotons");
      fThePhotonTreeDetected->Branch("EventID",     &fEventID,          "EventID/I");
      fThePhotonTreeDetected->Branch("Wavelength",  &fWavelength,       "Wavelength/F");
      fThePhotonTreeDetected->Branch("OpChannel",       &fOpChannel,            "OpChannel/I");
      fThePhotonTreeDetected->Branch("Time",        &fTime,             "Time/F");
    }
    
    if(fMakeOpDetsTree)
    {
      fTheOpDetTree    = tfs->make<TTree>("OpDets","OpDets");
      fTheOpDetTree->Branch("EventID",        &fEventID,          "EventID/I");
      fTheOpDetTree->Branch("OpChannel",          &fOpChannel,            "OpChannel/I");
      fTheOpDetTree->Branch("CountAll",       &fCountOpDetAll,      "CountAll/I");
      fTheOpDetTree->Branch("CountDetected",  &fCountOpDetDetected, "CountDetected/I");
      if(pvs->StoreReflected())
        fTheOpDetTree->Branch("CountReflDetected",  &fCountOpDetReflDetected, "CountReflDetected/I");
      fTheOpDetTree->Branch("Time",                   &fTime,                          "Time/F");
    }
    
    if(fMakeOpDetEventsTree)
    {
      fTheEventTree  = tfs->make<TTree>("OpDetEvents","OpDetEvents");
      fTheEventTree->Branch("EventID",      &fEventID,            "EventID/I");
      fTheEventTree->Branch("CountAll",     &fCountEventAll,     "CountAll/I");
      fTheEventTree->Branch("CountDetected",&fCountEventDetected,"CountDetected/I");
      if(pvs->StoreReflected())
        fTheOpDetTree->Branch("CountReflDetected",  &fCountOpDetReflDetected, "CountReflDetected/I");
      
    }
    
    //generating the tree for the light analysis:                                                                    
    if(fMakeLightAnalysisTree)
    {
      fLightAnalysisTree = tfs->make<TTree>("LightAnalysis","LightAnalysis");
      fLightAnalysisTree->Branch("RunNumber",&fRun);
      fLightAnalysisTree->Branch("EventID",&fEventID);
      fLightAnalysisTree->Branch("TrackID",&fTrackID);
      fLightAnalysisTree->Branch("PdgCode",&fpdg);
      fLightAnalysisTree->Branch("MotherTrackID",&fmotherTrackID);
      fLightAnalysisTree->Branch("Energy",&fEnergy);
      fLightAnalysisTree->Branch("dEdx",&fdEdx);
      fLightAnalysisTree->Branch("StepPrePositions",&fstepPrePositions);
      fLightAnalysisTree->Branch("StepPostPositions",&fstepPostPositions);
      fLightAnalysisTree->Branch("StepPreTimes",&fstepPreTimes);
      fLightAnalysisTree->Branch("StepPostTimes",&fstepPostTimes);
      fLightAnalysisTree->Branch("SignalsVUV",&fSignalsvuv);
      fLightAnalysisTree->Branch("SignalsVisible",&fSignalsvis);
      fLightAnalysisTree->Branch("Process",&fProcess);
      
      const int maxNtracks = 1000;
      fSignals_vuv.clear();
      fSignals_vuv.resize(maxNtracks);
      fSignals_vis.clear();
      fSignals_vis.resize(maxNtracks);
      for(size_t itrack=0; itrack!=maxNtracks; itrack++) {
        fSignals_vuv[itrack].resize(geo->NOpChannels());
        fSignals_vis[itrack].resize(geo->NOpChannels());
      }
      
      fstepPrePositions.clear();
      fstepPostPositions.clear();
      fstepPreTimes.clear();
      fstepPostTimes.clear();
      
    }
    
  }
  
  
  SimPhotonCounter::~SimPhotonCounter() 
  {
  }
  
  void SimPhotonCounter::endJob()
  {
    art::ServiceHandle<phot::PhotonVisibilityService> vis;
    
    if(vis->IsBuildJob())
    {
      vis->StoreLibrary();
    }
  }
  
  void SimPhotonCounter::analyze(art::Event const& evt)
  {
    
    // Lookup event ID from event
    art::EventNumber_t event = evt.id().event();
    fEventID=Int_t(event);
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;
    
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();
    
    // Service for determining opdet responses
    art::ServiceHandle<opdet::OpDetResponseInterface> odresponse;
    
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    
    //-------------------------stimation of dedx per trackID------------------------      
    
    //get the list of particles from this event       
    double totalEnergy_track[1000] = {0.};      
    if(fMakeLightAnalysisTree){
      const sim::ParticleList* plist = pi_serv? &(pi_serv->ParticleList()): nullptr;
      
      // loop over all sim::SimChannels in the event and make sure there are no             
      // sim::IDEs with trackID values that are not in the sim::ParticleList                
      std::vector<const sim::SimChannel*> sccol;
      //evt.getView(fG4ModuleLabel, sccol);                                                 
      evt.getView("largeant", sccol);
      double totalCharge=0.0;
      double totalEnergy=0.0;
      //loop over the sim channels collection  
      for(size_t sc = 0; sc < sccol.size(); ++sc){
        double numIDEs=0.0;
        double scCharge=0.0;
        double scEnergy=0.0;
        const auto & tdcidemap = sccol[sc]->TDCIDEMap();
        //loop over all of the tdc IDE map objects
        for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
          const std::vector<sim::IDE> idevec = (*mapitr).second;
          numIDEs += idevec.size();
          //go over all of the IDEs in a given simchannel
          for(size_t iv = 0; iv < idevec.size(); ++iv){
            if (plist) {
              if(plist->find( idevec[iv].trackID ) == plist->end()
                  && idevec[iv].trackID != sim::NoParticleId)
              {
                mf::LogWarning("LArG4Ana") << idevec[iv].trackID << " is not in particle list";
              }
            }
            if(idevec[iv].trackID < 0) continue;
            totalCharge +=idevec[iv].numElectrons;
            scCharge += idevec[iv].numElectrons;
            totalEnergy +=idevec[iv].energy;
            scEnergy += idevec[iv].energy;
            
            totalEnergy_track[idevec[iv].trackID] += idevec[iv].energy/3.;
          }
        }
      }
    }//End of if(fMakeLightAnalysisTree)

    
    if(!fUseLitePhotons)
    {
      //Get SimPhotonsCollection from Event
      sim::SimPhotonsCollection TheHitCollection = sim::SimListUtils::GetSimPhotonsCollection(evt,fInputModule);
      
      //Reset counters
      fCountEventAll=0;
      fCountEventDetected=0;
  
      if(TheHitCollection.size()>0)
      {
        if(fMakeLightAnalysisTree) {
          //resetting the signalt to save in the analysis tree per event
          const int maxNtracks = 1000;
          for(size_t itrack=0; itrack!=maxNtracks; itrack++) {
            for(size_t pmt_i=0; pmt_i!=geo->NOpChannels(); pmt_i++) {
              fSignals_vuv[itrack][pmt_i].clear();
              fSignals_vis[itrack][pmt_i].clear();
            }
          }
        }
        
        if(fVerbosity > 0) std::cout<<"Found OpDet hit collection of size "<< TheHitCollection.size()<<std::endl;
        if(TheHitCollection.size()>0)
        {
          for(sim::SimPhotonsCollection::const_iterator itOpDet=TheHitCollection.begin(); itOpDet!=TheHitCollection.end(); itOpDet++)
          {
            //Reset Counters
            fCountOpDetAll=0;
            fCountOpDetDetected=0;
            fCountOpDetReflDetected=0;
            //Reset t0 for visible light
            fT0_vis = 999.;
            
            //Get data from HitCollection entry
            fOpChannel=itOpDet->first;
            const sim::SimPhotons& TheHit=itOpDet->second;
            
            //std::cout<<"OpDet " << fOpChannel << " has size " << TheHit.size()<<std::endl;
            
            // Loop through OpDet phots.  
            //   Note we make the screen output decision outside the loop
            //   in order to avoid evaluating large numbers of unnecessary 
            //   if conditions. 
            
            if(fVerbosity > 3)
            {
              for(const sim::OnePhoton& Phot: TheHit)
              {
                // Calculate wavelength in nm
                fWavelength= odresponse->wavelength(Phot.Energy);
                
                //Get arrival time from phot
                fTime= Phot.Time;
                
                // Increment per OpDet counters and fill per phot trees
                fCountOpDetAll++;
                if(fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
                if(odresponse->detected(fOpChannel, Phot))
                {
                  if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
                  //only store direct direct light
                  if(!pvs->StoreReflected() || (pvs->StoreReflected() && fWavelength <200 )) 
                    fCountOpDetDetected++;
                  // reflected and shifted light is in visible range
                  else if(pvs->StoreReflected() && fWavelength >380 ) {
                    fCountOpDetReflDetected++;
                    // find the first visible arrival time 
                    if(pvs->StoreReflT0() && fTime < fT0_vis)
                      fT0_vis = fTime;  
                  }
                  
                  std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEventID<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 1 "<<std::endl;
                }
                else
                  std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEventID<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 0 "<<std::endl;
              }
            }
            else
            {
              for(const sim::OnePhoton& Phot: TheHit)
              {
                // Calculate wavelength in nm
                fWavelength= odresponse->wavelength(Phot.Energy);
                fTime= Phot.Time;    
                
                if(fMakeLightAnalysisTree) {
                  if(fWavelength <200)
                    fSignals_vuv[Phot.MotherTrackID][fOpChannel].push_back(fTime);
                  else
                    fSignals_vis[Phot.MotherTrackID][fOpChannel].push_back(fTime);
                  
                  initialPhotonPosition = Phot.InitialPosition;
                  finalPhotonPosition = Phot.FinalLocalPosition;
                }
                
                // Increment per OpDet counters and fill per phot trees
                fCountOpDetAll++;
                if(fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
                if(odresponse->detected(fOpChannel, Phot))
                {
                  if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
                  //only store direct/UV light
                  if(!pvs->StoreReflected() || (pvs->StoreReflected() && fWavelength <200 ))
                    fCountOpDetDetected++;
                  //shifted light is in visible range 
                  else if(pvs->StoreReflected() && fWavelength >380 ) { 
                    fCountOpDetReflDetected++;
                    // find the first visible arrival time
                    if(pvs->StoreReflT0() && fTime < fT0_vis )
                      fT0_vis= fTime;
                  }
                }
              }
            }
            
            // If this is a library building job, fill relevant entry
            art::ServiceHandle<phot::PhotonVisibilityService> pvs;
            if(pvs->IsBuildJob())
            {
              int VoxID; double NProd;
              pvs->RetrieveLightProd(VoxID, NProd);
              pvs->SetLibraryEntry(VoxID, fOpChannel, double(fCountOpDetDetected)/NProd);
              //store reflected light
              if(pvs->StoreReflected())
                pvs->SetLibraryEntry(VoxID, fOpChannel, double(fCountOpDetReflDetected)/NProd,true);
              //store reflected first arrival time 
              if(pvs->StoreReflT0())
                pvs->SetLibraryReflT0Entry(VoxID, fOpChannel, fT0_vis);
            }
            
            // Incremenent per event and fill Per OpDet trees      
            if(fMakeOpDetsTree) fTheOpDetTree->Fill();
            fCountEventAll+=fCountOpDetAll;
            fCountEventDetected+=fCountOpDetDetected;
            
            // Give per OpDet output
            if(fVerbosity >2) std::cout<<"OpDetResponseInterface PerOpDet : Event "<<fEventID<<" OpDet " << fOpChannel << " All " << fCountOpDetAll << " Det " <<fCountOpDetDetected<<std::endl; 
          }
          
          // Fill per event tree
          if(fMakeOpDetEventsTree) fTheEventTree->Fill();
          
          // Give per event output
          if(fVerbosity >1) std::cout<<"OpDetResponseInterface PerEvent : Event "<<fEventID<<" All " << fCountOpDetAll << " Det " <<fCountOpDetDetected<<std::endl;   
          
        }
        else
        {
          // if empty OpDet hit collection, 
          // add an empty record to the per event tree 
          if(fMakeOpDetEventsTree) fTheEventTree->Fill();
        }
        if(fMakeLightAnalysisTree) {
          std::cout<<"Building the analysis tree"<<std::endl;
          //---------------Building the analysis tree-----------:
          fRun = evt.run();
          std::vector<double> thisPrexyz;
          std::vector<double> thisPostxyz;
          
          art::Handle< std::vector<sim::MCTrack> > mctrackHandle;
          evt.getByLabel("mcreco",mctrackHandle);
          std::vector<sim::MCTrack> const& mctrackVec(*mctrackHandle);
          
          //loop over the particles (so over the tracks)                      
          fEnergy = mctrackVec[0][0].E(); 
          for(size_t i_p=0; i_p < mctrackVec.size(); i_p++){
            //resetting the vectors                                                                         
            fstepPrePositions.clear();
            fstepPostPositions.clear();
            fstepPreTimes.clear();
            fstepPostTimes.clear();
            fSignalsvuv.clear();
            fSignalsvis.clear();
            fdEdx = -1.;
            //filling the tree fields  
            fTrackID = mctrackVec[i_p].TrackID();
            fpdg = mctrackVec[i_p].PdgCode();
            fmotherTrackID = mctrackVec[i_p].MotherTrackID();
            fdEdx = totalEnergy_track[fTrackID];
            fSignalsvuv = fSignals_vuv[fTrackID];
            fSignalsvis = fSignals_vis[fTrackID];
            fProcess = mctrackVec[i_p].Process();
            //filling the center positions of each step                                   
            for(size_t i_s=1; i_s < mctrackVec[i_p].size(); i_s++){
              TVector3 const& vec1 = mctrackVec[i_p][i_s-1].Position().Vect();
              TVector3 const& vec2 = mctrackVec[i_p][i_s].Position().Vect();                                   
              thisPrexyz.clear();
              thisPrexyz.resize(3);
              thisPrexyz[0] = vec1.X();
              thisPrexyz[1] = vec1.Y();
              thisPrexyz[2] = vec1.Z();
              fstepPrePositions.push_back(thisPrexyz);
              thisPostxyz.clear();
              thisPostxyz.resize(3);
              thisPostxyz[0] = vec2.X();
              thisPostxyz[1] = vec2.Y();
              thisPostxyz[2] = vec2.Z();
              fstepPostPositions.push_back(thisPostxyz);
              fstepPreTimes.push_back(mctrackVec[i_p][i_s-1].T());
              fstepPostTimes.push_back(mctrackVec[i_p][i_s].T());
              //double stepL = (vec2-vec1).Mag(); 
              //std::cout<<"step length: "<<stepL<<std::endl;
              
            }
            //filling the tree per track
            fLightAnalysisTree->Fill();
          }  
        }
      }
    }

      
    if (fUseLitePhotons)
    {
      //Get SimPhotonsLite from Event
      art::Handle< std::vector<sim::SimPhotonsLite> > photonHandle; 
      evt.getByLabel("largeant", photonHandle);
      
      
      //Reset counters
      fCountEventAll=0;
      fCountEventDetected=0;
      
      if(fVerbosity > 0) std::cout<<"Found OpDet hit collection of size "<< (*photonHandle).size()<<std::endl;
      
      
      if((*photonHandle).size()>0)
      {
        
        for ( auto const& photon : (*photonHandle) )
        {
          //Get data from HitCollection entry
          fOpChannel=photon.OpChannel;
          std::map<int, int> PhotonsMap = photon.DetectedPhotons;
          
          //Reset Counters
          fCountOpDetAll=0;
          fCountOpDetDetected=0;
          
          if(fVerbosity > 3)
          {
            for(auto it = PhotonsMap.begin(); it!= PhotonsMap.end(); it++)
            {
              // Calculate wavelength in nm
              fWavelength= 128;
              
              //Get arrival time from phot
              fTime= it->first;
              //std::cout<<"Arrival time: " << fTime<<std::endl;
              
              for(int i = 0; i < it->second ; i++)
              {
                // Increment per OpDet counters and fill per phot trees
                fCountOpDetAll++;
                if(fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
                if(odresponse->detectedLite(fOpChannel))
                {
                  if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
                  fCountOpDetDetected++;
                  std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEventID<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 1 "<<std::endl;
                }
                else
                  std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEventID<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 0 "<<std::endl;
              }
            }
          }
          else
          {
            for(auto it = PhotonsMap.begin(); it!= PhotonsMap.end(); it++)
            {
              // Calculate wavelength in nm
              fWavelength= 128;
              fTime= it->first;    
              
              for(int i = 0; i < it->second; i++)
              {
                // Increment per OpDet counters and fill per phot trees
                fCountOpDetAll++;
                if(fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
                if(odresponse->detectedLite(fOpChannel))
                {
                  if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
                  fCountOpDetDetected++;
                }
              }
            }
          }
          
          // Incremenent per event and fill Per OpDet trees      
          if(fMakeOpDetsTree) fTheOpDetTree->Fill();
          fCountEventAll+=fCountOpDetAll;
          fCountEventDetected+=fCountOpDetDetected;
          
          // Give per OpDet output
          if(fVerbosity >2) std::cout<<"OpDetResponseInterface PerOpDet : Event "<<fEventID<<" OpDet " << fOpChannel << " All " << fCountOpDetAll << " Det " <<fCountOpDetDetected<<std::endl; 
        }
        // Fill per event tree
        if(fMakeOpDetEventsTree) fTheEventTree->Fill();
        
        // Give per event output
        if(fVerbosity >1) std::cout<<"OpDetResponseInterface PerEvent : Event "<<fEventID<<" All " << fCountOpDetAll << " Det " <<fCountOpDetDetected<<std::endl;   
        
      }
      else
      {
        // if empty OpDet hit collection, 
        // add an empty record to the per event tree 
        if(fMakeOpDetEventsTree) fTheEventTree->Fill();
      } 
    }
  }
}
namespace opdet{
  
  DEFINE_ART_MODULE(SimPhotonCounter)
    
}//end namespace opdet

