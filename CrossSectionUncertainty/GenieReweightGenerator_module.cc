#ifndef LAR1ND_NUANA
#define LAR1ND_NUANA value

/// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// Larsoft includes
#include "Geometry/Geometry.h"

// nusoft includes
#include "SimulationBase/MCFlux.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/GTruth.h"
#include "SimulationBase/MCParticle.h"

// root includes
#include <TTree.h>
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TROOT.h"

// c++ includes
#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>

// LArSoft includes                                                                                         
#include "alg/GenieReweightGeneratorAlg.h"




namespace rwgt{
  /// A module to extract the info from the Monte Carlo generator GENIE
  class GenieReweightGenerator : public art::EDAnalyzer {
  public:
    explicit GenieReweightGenerator(fhicl::ParameterSet const& pset);
    virtual ~GenieReweightGenerator(){};

    void beginJob();    
    // void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt);
    // void beginSubRun (const art::SubRun& subrun);
    void reset();

  private:

    // This alg does all of the computing, while this module interfaces with
    // art and packs all of the variables into the ttree
    GenieReweightGeneratorAlg fGenieReweightGeneratorAlg;



    //Variables needed from the .fcl file to get things out of the event
    std::string fGenieModuleLabel;

    std::vector<std::string> fWeights;
    std::vector<float> fWeightRangeSigma;
    unsigned int fRandSeed;
    int fNWeights;

    std::vector< std::vector<float> > reweightingSigmas;

    // #---------------------------------------------------------------
    // #This is the list of analysis variables needed for the ttree:
    // #---------------------------------------------------------------
    // These are the ttrees themselves:
    TTree* fTreeTot;    //This tree stores all the important information from an event

    std::vector<std::vector<float> > eventWeights;  

    
    // Genie Reweight vectors:
    std::vector< std::vector< float > > genieReweights;
    
    // #---------------------------------------------------------------
    // # End of the list of variables for the tree
    // #---------------------------------------------------------------
    

  }; // end of class GenieReweightGenerator
// }

  
  GenieReweightGenerator::GenieReweightGenerator(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fGenieModuleLabel (pset.get< std::string >              ("GenieModuleLabel"))
    , fWeights          (pset.get< std::vector<std::string> > ("Weights"))
    , fRandSeed         (pset.get< unsigned int >             ("RandSeed"))
    , fNWeights         (pset.get< int >                      ("NWeights"))
  {
  }


  void GenieReweightGenerator::beginJob(){

    reset();
    gROOT->ProcessLine(".L loadDictionaries.C+");

    // This function sets up the ttrees
    // get access to the TFile service  
    art::ServiceHandle<art::TFileService> tfs;

    fTreeTot = tfs->make<TTree>("EventsTot", "Event info for ALL types");


    fTreeTot->Branch("MultiWeight","MultiWeight",&eventWeights,32000,0);
    

    std::vector<reweight> reweights;
    fGenieReweightGeneratorAlg.parseWeights(fWeights, reweights);
    fGenieReweightGeneratorAlg.prepareSigmas(fNWeights, fRandSeed, reweightingSigmas);
    fGenieReweightGeneratorAlg.configureReWeight(reweights, reweightingSigmas);

    return;
  }

  void GenieReweightGenerator::reset(){
    
    // This method takes ANYTHING that goes into the ntuple and sets it to default.
    
    // Start by making sure the appropriate vectors are cleared and emptied:
    eventWeights.clear();
   
    return;
  }

  void GenieReweightGenerator::analyze(const art::Event& evt){

    this -> reset();
    
    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCTruth> > mclistGENIE;  
    art::Handle< std::vector<simb::MCFlux> > mcflux;
    art::Handle< std::vector<simb::GTruth> > mcgtruth;


    //actually go and get the stuff
    evt.getByLabel(fGenieModuleLabel,mclistGENIE);
    evt.getByLabel(fGenieModuleLabel,mcflux);
    evt.getByLabel(fGenieModuleLabel,mcgtruth);

    // contains the mctruth object from genie
    art::Ptr<simb::MCTruth> mc(mclistGENIE,0);
    
    // contains the mcflux object
    art::Ptr<simb::MCFlux > flux(mcflux,0);

    // contains the gtruth object
    art::Ptr<simb::GTruth > gtruth(mcgtruth,0);



    // Contains the neutrino info
    simb::MCNeutrino neutrino = mc -> GetNeutrino();
    
    std::cout << "The interaction info is: \n" 
              << "  gtruth->ftgtPDG................." << gtruth->ftgtPDG << "\n"
              << "  gtruth->ftgtZ..................." << gtruth->ftgtZ << "\n"
              << "  gtruth->ftgtA..................." << gtruth->ftgtA << "\n"
              << "  gtruth->fGint..................." << gtruth->fGint << "\n"
              << "  gtruth->fGscatter..............." << gtruth->fGscatter << "\n"
              << "  gtruth->fweight................." << gtruth->fweight << "\n"
              << "  neutrino.Mode()................." << neutrino.Mode() << "\n"
              << "  neutrino.InteractionType()......" << neutrino.InteractionType() << "\n"
              << "  neutrino.CCNC()................." << neutrino.CCNC() << "\n"
              << "  neutrino.Target()..............." << neutrino.Target() << "\n"
              << std::endl;

    // Now start packing up the variables to fill the tree
    // In general, have the algorithm do this:
    


    
    // Find a new weight for this event:
    // std::vector<std::vector<float> > weights;
    fGenieReweightGeneratorAlg.calcWeight(mc, gtruth,eventWeights);
    std::cout << "Printing weights for this event:\n";
    for (auto s : fWeights) std::cout << s << "\t";
    std::cout << "Total\n";
    for (unsigned int row = 0; row < eventWeights.front().size(); row ++){
    for (unsigned int column = 0; column < eventWeights.size(); column ++){
      std::cout << eventWeights[column][row] << "\t";
    }
    std::cout << std::endl;
    }


    // Fill the ttree:
    fTreeTot->Fill();


    return;
  }

  DEFINE_ART_MODULE(GenieReweightGenerator)


}

#endif