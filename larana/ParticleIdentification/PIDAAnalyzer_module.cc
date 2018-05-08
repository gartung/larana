////////////////////////////////////////////////////////////////////////
// Class:       PIDAAnalyzer
// Module Type: analyzer
// File:        PIDAAnalyzer_module.cc
//
// Generated at Sat Nov  1 23:10:27 2014 by Wesley Ketchum using artmod
// from cetpkgsupport v1_07_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "PIDAAlg.h"

namespace pid{
  class PIDAAnalyzer;
}

class pid::PIDAAnalyzer : public art::EDAnalyzer {
public:
  explicit PIDAAnalyzer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PIDAAnalyzer(PIDAAnalyzer const &) = delete;
  PIDAAnalyzer(PIDAAnalyzer &&) = delete;
  PIDAAnalyzer & operator = (PIDAAnalyzer const &) = delete;
  PIDAAnalyzer & operator = (PIDAAnalyzer &&) = delete;

  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p) ;
  void beginJob() override;

private:

  std::string fCaloModuleLabel;
  PIDAAlg     fPIDAAlg;

};


pid::PIDAAnalyzer::PIDAAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCaloModuleLabel(p.get<std::string>("CaloModuleLabel")),
  fPIDAAlg(p.get<fhicl::ParameterSet>("PIDAAlg"))
{}

void pid::PIDAAnalyzer::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;

  std::vector<TH1F*> kde_hists;
  for(size_t i_b=0; i_b<fPIDAAlg.getNKDEBandwidths(); i_b++){
    std::stringstream hname;
    hname << "hkde_" << i_b;
    kde_hists.push_back(tfs->make<TH1F>(hname.str().c_str(),"PIDA KDE Distribution",100,0,30));
  }

  fPIDAAlg.SetPIDATree(tfs->make<TTree>("pida","PIDAPropertiesTree"),
		       tfs->make<TH1F>("hvals","PIDA Distribution",100,0,30),
		       kde_hists);
		      
}

void pid::PIDAAnalyzer::analyze(art::Event const & e)
{
  art::Handle< std::vector<anab::Calorimetry> > caloHandle;
  e.getByLabel(fCaloModuleLabel,caloHandle);
  std::vector<anab::Calorimetry> const& caloVector(*caloHandle);

  for(size_t i_calo=0; i_calo<caloVector.size(); i_calo++){
    fPIDAAlg.FillPIDATree(e.run(),e.event(),i_calo,caloVector[i_calo]);
  }

}

void pid::PIDAAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  fCaloModuleLabel = p.get<std::string>("CaloModuleLabel");
  fPIDAAlg.reconfigure(p.get<fhicl::ParameterSet>("PIDAAlg"));
}

DEFINE_ART_MODULE(pid::PIDAAnalyzer)
