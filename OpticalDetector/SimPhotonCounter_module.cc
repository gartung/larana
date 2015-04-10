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
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "OpticalDetector/OpDetResponseInterface.h"
#include "Simulation/SimListUtils.h"
#include "Simulation/sim.h"
#include "Simulation/LArG4Parameters.h"

// ROOT includes
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>

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
     
	//Time ahack
	  void clearVectors(int OpCh, int Voxels);

    private:
      
      // Trees to output

      TTree * fThePhotonTreeAll;
      TTree * fThePhotonTreeDetected;
      TTree * fTheOpDetTree;
      TTree * fTheEventTree;

	// Timing Tree --ahack
	  TTree * fTheTimeTree;
	  //TH1F * fTheTimeHist ; 
 	//  TObjArray *timeHists ;
	 // TObjArray *distHists ;

     TH1F *timeHist;//=new TH1F(histname,"",200,0,80);
     TH1F *distHist; //=new TH1F(disthistname,"",200,0,80);
	  
	  
      // Parameters to read in

      std::string fInputModule;      // Input tag for OpDet collection

      int fVerbosity;                // Level of output to write to std::out

      bool fMakeDetectedPhotonsTree; //
      bool fMakeAllPhotonsTree;      //
      bool fMakeOpDetsTree;         // Switches to turn on or off each output
      bool fMakeOpDetEventsTree;          //
	  bool fMakeTimeTree = false; 
      
      float fQE;                     // Quantum efficiency of tube

      float fWavelengthCutLow;       // Sensitive wavelength range 
      float fWavelengthCutHigh;      // 


      

      // Data to store in trees

      Float_t fWavelength;
      Float_t fTime;
      Int_t fCount;
      Int_t fCountOpDetAll;
      Int_t fCountOpDetDetected;

      Int_t fCountEventAll;
      Int_t fCountEventDetected;
      
      Int_t fEventID;
      Int_t fOpChannel;


	//Time Tree --ahack
	  std::vector<std::vector<std::vector<float>>> fLookupTime ;
	  std::vector<std::vector<std::vector<float>>> fLookupDist ;
	  Int_t	fVoxel;
	  Float_t fDist;

	//Corey's fix
	  std::map<int,int> voxelToIndex ;
	  std::vector<int> fVoxelList ;

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

    //fQE=                       pset.get<double>("QuantumEfficiency");
    //fWavelengthCutLow=         pset.get<double>("WavelengthCutLow");
    //fWavelengthCutHigh=        pset.get<double>("WavelengthCutHigh");
    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);    
  }
  

  void SimPhotonCounter::beginJob()
  {
    // Get file service to store trees
    art::ServiceHandle<art::TFileService> tfs;

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
	fTheOpDetTree->Branch("Time",  			&fTime,				 "Time/F");
      }
    
    if(fMakeOpDetEventsTree)
      {
	fTheEventTree  = tfs->make<TTree>("OpDetEvents","OpDetEvents");
	fTheEventTree->Branch("EventID",      &fEventID,            "EventID/I");
	fTheEventTree->Branch("CountAll",     &fCountEventAll,     "CountAll/I");
	fTheEventTree->Branch("CountDetected",&fCountEventDetected,"CountDetected/I");
      }

	//ahack
     if(fMakeTimeTree)
       {
     fTheTimeTree = tfs->make<TTree>("PerVoxOpCh","PerVoxOpCh");
	 fTheTimeTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
	 fTheTimeTree->Branch("Voxel", &fVoxel, "Voxel/I");

//	 fTheTimeHist = tfs->make<TH1F>("TimeHist","Time per vox per opch",100,0,120) ;
	// timeHists = tfs->make<TObjArray>();
	// distHists = tfs->make<TObjArray>();
       timeHist=new TH1F("timehist","",200,0,80);
       distHist=new TH1F("disthist","",200,0,80);
	  
	  
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

	//Time --ahack
	void SimPhotonCounter::clearVectors(int OpCh, int Voxels){

		fLookupTime.clear();
		fLookupTime.resize(OpCh); 	

		fLookupDist.clear();
		fLookupDist.resize(OpCh); 	

	 	for(int i=0; i < OpCh; i++){
    	   fLookupTime[i].resize(Voxels);  
    	   fLookupDist[i].resize(Voxels);  
			for(int j=0; j < Voxels; j++){
			  fLookupTime[i][j].resize(0) ;
			  fLookupDist[i][j].resize(0) ;
				}
         }
		}


  void SimPhotonCounter::analyze(art::Event const& evt)
  {

    // Setup random number generator (for QE sampling)
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);

    // Lookup event ID from event
    art::EventNumber_t event = evt.id().event();
    fEventID=Int_t(event);

 	//Time --ahack
	art::ServiceHandle<phot::PhotonVisibilityService> pvs;



    art::ServiceHandle<sim::LArG4Parameters> lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();


    // Service for determining opdet responses
    art::ServiceHandle<opdet::OpDetResponseInterface> odresponse;


    if(!fUseLitePhotons)
    {
    //Get SimPhotonsCollection from Event
    sim::SimPhotonsCollection TheHitCollection = sim::SimListUtils::GetSimPhotonsCollection(evt,fInputModule);

    //Reset counters
    fCountEventAll=0;
    fCountEventDetected=0;

	//ahack
	auto VoxelDef = pvs->GetVoxelDef();
	auto Voxels = VoxelDef.GetNVoxels();

       //std::cout << " voxels: " << Voxels << " Def " << VoxelDef << std::endl;
       std::cout << " voxels: " << Voxels << " Def "  << std::endl;
   // auto Voxels = fVoxelList.size();
  //  clearVectors(TheHitCollection.size(),Voxels);
	//std::cout<<"Voxels! "<<Voxels<<std::endl;
	
//ahack
//	TH1F *h[TheHitCollection.size()][Voxels]; 
//	for(Int_t o=0; o< TheHitCollection.size(); o++){
//		  h[o][fEventID-1] = new TH1F(Form("h%f_%f",o,fEventID-1),"Title",100,0,120);
//	}

    if(fVerbosity > 0) std::cout<<"Found OpDet hit collection of size "<< TheHitCollection.size()<<std::endl;
    if(TheHitCollection.size()>0)
      {
	for(sim::SimPhotonsCollection::const_iterator itOpDet=TheHitCollection.begin(); itOpDet!=TheHitCollection.end(); itOpDet++)
	  {
	    //Reset Counters
	    fCountOpDetAll=0;
	    fCountOpDetDetected=0;

	    //Get data from HitCollection entry
	    fOpChannel=itOpDet->first;
	    const sim::SimPhotons& TheHit=itOpDet->second;
	     
//	        std::cout<<"OpDet " << fOpChannel << " has size " << TheHit.size()<<std::endl;
	    
	    // Loop through OpDet phots.  
	    //   Note we make the screen output decision outside the loop
	    //   in order to avoid evaluating large numbers of unnecessary 
	    //   if conditions. 
   
	if(fMakeTimeTree)
	  {	
	    //ahack ->Fill Histogram
	   char *histname = new char[10];
	   char *disthistname = new char[10];	 
	
	    sprintf(histname, "hTime_%i_%i",fOpChannel,fEventID-1);
	    timeHist->Reset();
	    timeHist->SetName(histname);

	    sprintf(disthistname, "hDist_%i_%i",fOpChannel,fEventID-1);
	    distHist->Reset();
            distHist->SetName(disthistname); 
	   
	  } 
		
	    if(fVerbosity > 3)
	      {
		for(const sim::OnePhoton& Phot: TheHit)
		  {
		    // Calculate wavelength in nm
                    fWavelength= odresponse->wavelength(Phot.Energy);

		    //Get arrival time from phot
		    fTime = Phot.Time;

		    std::cout<<"Arrival time: " << fTime<<std::endl;
		    
		    // Increment per OpDet counters and fill per phot trees
		    fCountOpDetAll++;
		    if(fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
                    if(odresponse->detected(fOpChannel, Phot))
		      {
			if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
			fCountOpDetDetected++;
			std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEventID<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 1 "<<std::endl;
		      }
		    else
		      std::cout<<"OpDetResponseInterface PerPhoton : Event "<<fEventID<<" OpChannel " <<fOpChannel << " Wavelength " << fWavelength << " Detected 0 "<<std::endl;
		  }
	      }
	    else
	      {
		//ahack
		// having segfault issues with this;not sure why. Using 1/c below.
		//	quantity is 0.0333564 
		//	float speedOfLight = 29.9792 ; //in cm/ns
	//	std::cout<<"Voxel, OpChannel: "<<fEventID -1<<", "<<fOpChannel<<std::endl;
		for(const sim::OnePhoton& Phot: TheHit)
		  {
		    // Calculate wavelength in nm
                    fWavelength= odresponse->wavelength(Phot.Energy);
		    fTime= Phot.Time;		

		 
		    // Increment per OpDet counters and fill per phot trees
		    fCountOpDetAll++;
		    if(fMakeAllPhotonsTree) fThePhotonTreeAll->Fill();
                    if(odresponse->detected(fOpChannel, Phot))
		      {
			if(fMakeDetectedPhotonsTree) fThePhotonTreeDetected->Fill();
			fCountOpDetDetected++;
		      }

//			std::cout<<"Stopping before filling time and dist vecs"<<std::endl;
			//Fill time&dist assoc with an OpChannel and voxel --ahack
			if(pvs->IsBuildJob()) //&& voxelToIndex.find(fEventID-1)!= voxelToIndex.end() )
                   {
//				std::cout<<"Adding to the vector of vector...etc.  "<<std::endl;
		if(fMakeTimeTree)
		  {	
		      //Get distance of photon's start point from point located on OpDet and convert
		   //mm to cm (LArG4)
		   fDist =0.0333564*0.1*pow(pow(Phot.InitialPosition.X() - Phot.FinalPosition.X(),2)
			         	  + pow(Phot.InitialPosition.Y() - Phot.FinalPosition.Y(),2) 
					  + pow(Phot.InitialPosition.Z() - Phot.FinalPosition.Z(),2),0.5);  

	       //		std::cout<<"Initial Pos :"<<Phot.InitialPosition.X() <<std::endl;
	       //				 <<" "<<Phot.InitialPosition.Y()<<" "<<Phot.InitialPosition.Z()
	       //		std::cout<<"\nFinal Pos: "<<Phot.FinalPosition.X()/10
	       //				 <<" "<<Phot.FinalPosition.Y()/10<<" "<<Phot.FinalPosition.Z()/10<<std::endl;
		    timeHist->Fill(fTime);
		    distHist->Fill(fDist);
		    
		  } 
               //  fLookupTime[fOpChannel][fEventID-1].push_back(fTime) ;
               //  fLookupDist[fOpChannel][fEventID-1].push_back(fDist) ;

				 //std::cout<<"\nTime: "<<fTime<<", size while filling: ["<<fOpChannel<<"]["<<VoxID<<"] "<<fLookupTime[fOpChannel][VoxID].size();

					}

		  }
	      }
	  
	    if(fMakeTimeTree)
		    {
		    timeHist->Write();     // all objects from timeHists are written 
		    distHist->Write();
		    }
	  	      
	    // If this is a library building job, fill relevant entry
//	    art::ServiceHandle<phot::PhotonVisibilityService> pvs;
	    if(pvs->IsBuildJob())
	      {
		int VoxID; double NProd;
		pvs->RetrieveLightProd(VoxID, NProd);
		pvs->SetLibraryEntry(VoxID, fOpChannel, double(fCountOpDetDetected)/NProd);		
	      }


	    // Incremenent per event and fill Per OpDet trees	    
	    if(fMakeOpDetsTree) fTheOpDetTree->Fill();
	    fCountEventAll+=fCountOpDetAll;
	    fCountEventDetected+=fCountOpDetDetected;

	    // Give per OpDet output
	    if(fVerbosity >2) std::cout<<"OpDetResponseInterface PerOpDet : Event "<<fEventID<<" OpDet " << fOpChannel << " All " << fCountOpDetAll << " Det " <<fCountOpDetDetected<<std::endl; 
	  }

	
//      for(size_t o = 0; o < TheHitCollection.size(); o++){
// 	   //only want to store voxel for matching event ID because 1 event is generated in 1 voxel
// 	   int v = fEventID - 1 ;
// 
// //	   std::cout<<"\nThe timesize after filling ["<<o<<"]["<<v<<"]: "<<fLookupTime[o][v].size();
// //		std::cout<<"i hate larsoft "<<std::endl;
// 	   if(pvs->IsBuildJob()){// && voxelToIndex.find(fEventID-1) !=voxelToIndex.end() ){
// //	   	  int v = voxelToIndex[fEventID-1];
// 	   if(fLookupTime[int(o)][v].size() > 0. ){
// 	 	  sprintf(histname, "hTime_%i_%i",int(o),fEventID-1);
// 		  TH1F *timeHist=new TH1F(histname,"",200,0,80);
// 
// 	 	  sprintf(disthistname, "hDist_%i_%i",int(o),fEventID-1);
// 		  TH1F *distHist=new TH1F(disthistname,"",200,0,80);
// 
// 		  timeHists->AddLast(timeHist);
// 		  distHists->AddLast(distHist);
// 		  //LookupTime and LookupDist should have same number of elements (distance and time are saved
// 		  //per photon, per voxel, per opdet
// 		  for(size_t t=0; t< fLookupTime[int(o)][v].size(); t++){
// 			timeHist->Fill(fLookupTime[int(o)][v][t]);
// 			distHist->Fill(fLookupDist[int(o)][v][t]);
// 		 	 }
// //	  for(size_t t=0; t< fLookupTime[int(o)][voxelToIndex[v]].size(); t++){
// //            timeHist->Fill(fLookupTime[int(o)][voxelToIndex[v]][t]);
// //            distHist->Fill(fLookupDist[int(o)][voxelToIndex[v]][t]);
// //              }          
// 		  
// 	//	  for(size_t d=0; d< fLookupDist[int(o)][v].size(); d++)
// 		 	 
// 		  
// 
// 		   timeHist->Write();     // all objects from timeHists are written 
// 		   distHist->Write();
// 		 	} 
// 		 }//if IsBuildJob
// 		else{
// 			continue;
// 			std::cout<<"we should never get here..."<<std::endl;
// 		}
// 
//      }       


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
    else
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
		    fTime= it->first*2;
		    std::cout<<"Arrival time: " << fTime<<std::endl;
		   
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
		      fTime= it->first*2;		
   
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
