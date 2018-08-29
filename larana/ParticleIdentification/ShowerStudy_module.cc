////////////////////////////////////////////////////////////////////////
//
// file ShowerStudy_module.cc
//
// j.pillow@warwick.ac.uk
//
///////////////////////////////////////////////////////////////////////

// Generic C++ includes
#include <iostream>
#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h" // Shower energy finder
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"

#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TH2D.h"
#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<Double_t> >+;
#endif

namespace protoDUNE
{
  class ShowerStudyAna : public art::EDAnalyzer
  {
    public:
    
      explicit ShowerStudyAna (fhicl::ParameterSet const & pset);
    
      void reconfigure(fhicl::ParameterSet const & pset);
    
      void analyze ( art::Event const & evt) override;
    
      void beginJob();
    
      void PFenergyReco ( art::FindManyP<recob::Cluster>             const & PFtoClust,
                          art::FindManyP<recob::Cluster>             const & HitToClust,
                          std::vector< art::Ptr<recob::Hit> >        const & hitList,
                          std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                          double const & t0,
                          bool   const & cheat );
    
      void TrackShowerEnergyReco ( art::FindManyP<recob::Hit>             const & TSToHit,
                                   std::vector< art::Ptr<recob::Shower> > const & TSList,
                                   double const & t0,
                                   bool   const & emshower,
                                   bool   const & track,
                                   bool   const & cheat );
    
      void TrackShowerEnergyReco ( art::FindManyP<recob::Hit>             const & TSToHit,
                                   std::vector< art::Ptr<recob::Track> >  const & TSList,
                                   double const & t0,
                                   bool   const & emshower,
                                   bool   const & track,
                                   bool   const & cheat );
    
      void PFenergyRecoByFlow ( art::FindManyP<recob::Cluster>             const & PFtoClust,
                                art::FindManyP<recob::Cluster>             const & HitToClust,
                                std::vector< art::Ptr<recob::Hit> >        const & hitList,
                                std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                double const & t0,
                                bool   const & cheat );
    
      std::vector<double> DaughterEnergy ( std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                           art::FindManyP<recob::Cluster>             const & PFtoClust,
                                           std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > const & clustersToHits,
                                           std::vector<size_t> const & daughters,
                                           double const & t0);
    

    
    private:
    
      // Generic Variables
      double t0;
      double initialE;
      double frontE;
      double frontX;
      double frontY;
      double frontZ;
    
      double recombCorrection { 1.0/0.63 };
      double ADCtoGeV         { 23.6e-9 / 4.966e-3 };
    
      const simb::MCParticle* pPartPrimary;
    
      // Algorithms for Energy Reco
      shower::ShowerEnergyAlg fShowerEnergyAlg;
      calo::CalorimetryAlg    fCaloAlg;
    
      // Root Tree
      TTree* showerTree;
    
      // Root Tree Variables
    
      std::vector<double> MCIDE;
    
      std::vector<double> PFTotCheat;
      std::vector<double> PFShowersTotCheat;
      std::vector<double> PFTracksTotCheat;
   
      std::vector<double> PFTot;
      std::vector<double> PFShowersTot;
      std::vector<double> PFTracksTot;
    
      std::vector< std::vector<double> > PFShowersCheat;
      std::vector< std::vector<double> > PFTracksCheat;
    
      std::vector< std::vector<double> > PFShowers;
      std::vector< std::vector<double> > PFTracks;
    
      std::vector< std::vector<double> > emShower;
      std::vector< std::vector<double> > pandoraTrack;
      std::vector< std::vector<double> > pandoraShower;
      std::vector< std::vector<double> > pandoraShowerCheat;
    
      std::vector<double> emShowerTotE;
      std::vector<double> pandoraTrackTotE;
      std::vector<double> pandoraShowerTotE;
      std::vector<double> pandoraShowerCheatTotE;
    
      std::vector< std::vector<double> > EfromPFlow;
      std::vector< std::vector<double> > EfromPFlowCheat;
    
      std::vector< int > PFlowPDGPrimary;
      std::vector< int > PFlowPDGPrimaryCheat;
    
      std::vector< std::vector<double> > spacePoints;
      std::vector< std::vector<double> > spacePointsCheat;
    
      std::string MCModuleLabel;
      std::string MCPartModLabel;
      std::string emshowerLabel;
      std::string pandoraShowerLabel;
      std::string pandoraTrackLabel;
      std::string pandoraShowerLabelCheat;
      std::string HitLabel;
      std::string PFParticleLabel;
      std::string spacePointLabel;
      std::string spacePointLabelCheat;
      std::string clusterLabel;
      std::string PFParticleLabelCheat;
      std::string clusterLabelCheat;
    
  }; // class ShowerStudyAna : public art::EDAnalyzer
  
} // namespace protoDUNE

//---------------------------------------------------------------------------------------------------------

protoDUNE::ShowerStudyAna::ShowerStudyAna(fhicl::ParameterSet const & pset):
  EDAnalyzer(pset),
  fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
  fCaloAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  this->reconfigure(pset);
} // ShowerStudyAna::ShowerStudyAna(fhicl::ParameterSet const& pset)

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  
  showerTree = tfs->make<TTree>( "showerTree", "showerTree" );
  
  showerTree->Branch( "initialE", &initialE );
  showerTree->Branch( "frontE",   &frontE   );
  showerTree->Branch( "frontX",   &frontX   );
  showerTree->Branch( "frontY",   &frontY   );
  showerTree->Branch( "frontZ",   &frontZ   );
  
  showerTree->Branch( "MCIDE",   &MCIDE   );
  
  showerTree->Branch( "PFTotCheat",        &PFTotCheat        );
  showerTree->Branch( "PFShowersTotCheat", &PFShowersTotCheat );
  showerTree->Branch( "PFTracksTotCheat",  &PFTracksTotCheat  );
  
  showerTree->Branch( "PFTot",        &PFTot        );
  showerTree->Branch( "PFShowersTot", &PFShowersTot );
  showerTree->Branch( "PFTracksTot",  &PFTracksTot  );
  
  showerTree->Branch( "PFShowersCheat", &PFShowersCheat );
  showerTree->Branch( "PFTracksCheat",  &PFTracksCheat  );
  
  showerTree->Branch( "PFShowers", &PFShowers );
  showerTree->Branch( "PFTracks",  &PFTracks  );
  
  showerTree->Branch( "emShower",           &emShower           );
  showerTree->Branch( "pandoraTrack",       &pandoraTrack       );
  showerTree->Branch( "pandoraShower",      &pandoraShower      );
  showerTree->Branch( "pandoraShowerCheat", &pandoraShowerCheat );

  showerTree->Branch( "emShowerTotE",           &emShowerTotE           );
  showerTree->Branch( "pandoraTrackTotE",       &pandoraTrackTotE       );
  showerTree->Branch( "pandoraShowerTotE",      &pandoraShowerTotE      );
  showerTree->Branch( "pandoraShowerCheatTotE", &pandoraShowerCheatTotE );
  
  showerTree->Branch( "EfromPFlow",      &EfromPFlow      );
  showerTree->Branch( "EfromPFlowCheat", &EfromPFlowCheat );
  
  showerTree->Branch( "PFlowPDGPrimary",      &PFlowPDGPrimary      );
  showerTree->Branch( "PFlowPDGPrimaryCheat", &PFlowPDGPrimaryCheat );
  
  showerTree->Branch( "spacePoints",      &spacePoints      );
  showerTree->Branch( "spacePointsCheat", &spacePointsCheat );

}

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::reconfigure(fhicl::ParameterSet const & p)
{ // ==================================================================
  // Data labels
  // ==================================================================

  MCModuleLabel           = p.get< std::string >("MCModuleLabel");            // "generator";
  MCPartModLabel          = p.get< std::string >("MCPartModLabel");           // "largeant";
  emshowerLabel           = p.get< std::string >("emshowerLabel");            // "emshower";
  pandoraShowerLabel      = p.get< std::string >("pandoraShowerLabel");       // "pandoraShower";
  pandoraTrackLabel       = p.get< std::string >("pandoraTrackLabel");        // "pandoraTrack";
  pandoraShowerLabelCheat = p.get< std::string >("pandoraShowerLabelCheat");  // "cheatPandoraShower";
  HitLabel                = p.get< std::string >("HitLabel");                 // "linecluster";
  PFParticleLabel         = p.get< std::string >("PFParticleLabel");          // "pandora";
  spacePointLabel         = p.get< std::string >("spacePointLabel");          // "pandora";
  spacePointLabelCheat    = p.get< std::string >("spacePointLabelCheat");     // "cheatPandora";
  clusterLabel            = p.get< std::string >("clusterLabel");             // "pandora";
  PFParticleLabelCheat    = p.get< std::string >("PFParticleLabelCheat");     // "cheatPandora";
  clusterLabelCheat       = p.get< std::string >("clusterLabelCheat");        // "cheatPandora";

}

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::analyze( art::Event const & evt )
{
  

  // ==================================================================
  // Clear initialised vectors
  // ==================================================================
  
  MCIDE.clear();
  
  PFTotCheat.clear();
  PFShowersTotCheat.clear();
  PFTracksTotCheat.clear();

  PFTot.clear();
  PFShowersTot.clear();
  PFTracksTot.clear();

  PFShowersCheat.clear();
  PFTracksCheat.clear();

  PFShowers.clear();
  PFTracks.clear();
  
  emShower.clear();
  pandoraTrack.clear();
  pandoraShower.clear();
  pandoraShowerCheat.clear();
  
  emShowerTotE.clear();
  pandoraTrackTotE.clear();
  pandoraShowerTotE.clear();
  pandoraShowerCheatTotE.clear();
  
  EfromPFlow.clear();
  EfromPFlowCheat.clear();
  
  PFlowPDGPrimary.clear();
  PFlowPDGPrimaryCheat.clear();
  
  spacePoints.clear();
  spacePointsCheat.clear();
  
  // ==================================================================
  // Get the data and service handles
  // ==================================================================
  
  // MCTruth
  art::Handle< std::vector<simb::MCTruth> > mcHandle;
  std::vector< art::Ptr<simb::MCTruth> > mcList;
  if ( evt.getByLabel( MCModuleLabel, mcHandle ) ) art::fill_ptr_vector( mcList, mcHandle );
  
  // MCParticle
  art::Handle< std::vector<simb::MCParticle> > mcPartHandle;
  std::vector< art::Ptr<simb::MCParticle> > mcPartList;
  if ( evt.getByLabel( MCPartModLabel, mcPartHandle ) ) art::fill_ptr_vector( mcPartList, mcPartHandle );
  
  // Hits linecluster
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector< art::Ptr<recob::Hit> > hitList;
  if ( evt.getByLabel( HitLabel, hitHandle ) ) art::fill_ptr_vector( hitList, hitHandle );
  
  // PFParticles
  art::Handle< std::vector<recob::PFParticle> > PFPartHandle;
  std::vector< art::Ptr<recob::PFParticle> > PFPartList;
  if ( evt.getByLabel( PFParticleLabel, PFPartHandle ) ) art::fill_ptr_vector( PFPartList, PFPartHandle );
  
  // Clusters
  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  std::vector< art::Ptr<recob::Cluster> > clusterList;
  if ( evt.getByLabel( clusterLabel, clusterHandle ) ) art::fill_ptr_vector( clusterList, clusterHandle );
  
  // Cheated PFParticles
  art::Handle< std::vector<recob::PFParticle> > PFPartHandleCheat;
  std::vector< art::Ptr<recob::PFParticle> > PFPartListCheat;
  if ( evt.getByLabel( PFParticleLabelCheat, PFPartHandleCheat ) ) art::fill_ptr_vector( PFPartListCheat, PFPartHandleCheat );
  
  // Cheated Clusters
  art::Handle< std::vector<recob::Cluster> > clusterHandleCheat;
  std::vector< art::Ptr<recob::Cluster> > clusterListCheat;
  if ( evt.getByLabel( clusterLabelCheat, clusterHandleCheat ) ) art::fill_ptr_vector( clusterListCheat, clusterHandleCheat );
  
  // pandoraShower
  art::Handle< std::vector<recob::Shower> > pandoraHandle;
  std::vector< art::Ptr<recob::Shower> > pandoraList;
  if ( evt.getByLabel( pandoraShowerLabel, pandoraHandle ) ) art::fill_ptr_vector( pandoraList, pandoraHandle );
  
  // pandoraShower Cheated
  art::Handle< std::vector<recob::Shower> > pandoraHandleCheat;
  std::vector< art::Ptr<recob::Shower> > pandoraListCheat;
  if ( evt.getByLabel( pandoraShowerLabelCheat, pandoraHandleCheat ) ) art::fill_ptr_vector( pandoraListCheat, pandoraHandleCheat );
  
  // pandoraTrack
  art::Handle< std::vector<recob::Track> > pandoraTrackHandle;
  std::vector< art::Ptr<recob::Track> > pandoraTrackList;
  if ( evt.getByLabel( pandoraTrackLabel, pandoraTrackHandle ) ) art::fill_ptr_vector( pandoraTrackList, pandoraTrackHandle );
  
  // emshower
  art::Handle< std::vector<recob::Shower> > emshowerHandle;
  std::vector< art::Ptr<recob::Shower> > emshowerList;
  if ( evt.getByLabel( emshowerLabel, emshowerHandle ) ) art::fill_ptr_vector( emshowerList, emshowerHandle );
  
  // spacepoints
  art::Handle< std::vector<recob::SpacePoint> > spacePointHandle;
  std::vector< art::Ptr<recob::SpacePoint> > spacePointList;
  if ( evt.getByLabel( spacePointLabel, spacePointHandle ) ) art::fill_ptr_vector( spacePointList, spacePointHandle );
  
  // spacepoints Cheated
  art::Handle< std::vector<recob::SpacePoint> > spacePointHandleCheat;
  std::vector< art::Ptr<recob::SpacePoint> > spacePointListCheat;
  if ( evt.getByLabel( spacePointLabelCheat, spacePointHandleCheat ) ) art::fill_ptr_vector( spacePointListCheat, spacePointHandleCheat );
  
  
  // ==================================================================
  // Get the initial energy of the original MC Particle and T0
  // ==================================================================
  
  // Initial energy
  art::Ptr<simb::MCTruth> mcTruthProto = mcList[0];
  const simb::MCParticle& partp( mcTruthProto->GetParticle(0) );
  initialE = partp.E();

  // Get Trajectory shizzle
  auto traj = partp.Trajectory();

  // T0
  auto const *detprop = lar::providerFrom< detinfo::DetectorPropertiesService >();
  t0 = detprop->TriggerOffset();
//  double g4Init = 0;
  // Run through the geant tracks to get the tracjectory stuff. Match by initial energy
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
  sim::ParticleList::const_iterator itPart = plist.begin(), pend = plist.end();
  for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
    const simb::MCParticle* pPart = (itPart++)->second;
    if (!pPart) {
      throw art::Exception(art::errors::LogicError) << "GEANT particle #" << iPart << " returned a null pointer";
    }
    auto trajG = pPart->Trajectory();
//    std::cout << "Traj energy: " << trajG[0].second.E() << "\n"
//              << ".E() energy: " << pPart->E() << "\n";
    if ( fabs(trajG[0].second.E() - initialE) < 0.00001 ){
//      std::cout << "\ntrajG.size(): " << trajG.size() << " | Length: " << trajG.TotalLength() << " | Energy: " << trajG[0].second.E() << "\n";
//      g4Init = trajG[0].second.E();
      pPartPrimary = pPart;
      break;
    }
    
  }
  
  
  // Find the energy at the front of the detector, Z = -0.49375
  art::ServiceHandle<geo::Geometry> geom;
  double EFF = 0;
  double XFF = 0;
  double YFF = 0;
  double ZFF = 0;
  auto trajP = pPartPrimary->Trajectory();
  for ( auto it = trajP.begin() ; it != trajP.end() ; it++ ) {
    double point[3] = { it->first.X(), it->first.Y(), it->first.Z() };
    geo::TPCID idtpc = geom->FindTPCAtPosition(point);
    if ( geom->HasTPC(idtpc) ) {
      EFF = it->second.E();
      XFF = it->first.X();
      YFF = it->first.Y();
      ZFF = it->first.Z();
      std::cout << "X: " << it->first.X() << " | Y: " << it->first.Y() << " | Z: " << it->first.Z() << "\n";
      break;
    }
  }
  frontE = EFF;
  frontX = XFF;
  frontY = YFF;
  frontZ = ZFF;
//  std::cout << "G4 Initial Energy: " << g4Init << "\n"
//            << "MC Initial Energy: " << initialE << "\n"
//            << "Front Face Energy: " << EFF << "\n";
 
  // ==================================================================
  // Associations
  // ==================================================================

  art::FindManyP< recob::Hit > pandoraToHit      ( pandoraHandle,      evt, pandoraShowerLabel      );
  art::FindManyP< recob::Hit > pandoraToHitCheat ( pandoraHandleCheat, evt, pandoraShowerLabelCheat );
  art::FindManyP< recob::Hit > emshowerToHit     ( emshowerHandle,     evt, emshowerLabel           );
  
  art::FindManyP< recob::Hit > pandoraTrackToHit ( pandoraTrackHandle, evt, pandoraTrackLabel      );

  art::FindManyP< recob::Cluster > PFtoClust       ( PFPartHandle,      evt, PFParticleLabel      );
  art::FindManyP< recob::Cluster > HitToClust      ( hitHandle,         evt, PFParticleLabel      );
  art::FindManyP< recob::Shower  > PFtoShower      ( PFPartHandle,      evt, PFParticleLabel      );
  
  art::FindManyP< recob::Cluster > PFtoClustCheat  ( PFPartHandleCheat, evt, PFParticleLabelCheat );
  art::FindManyP< recob::Cluster > HitToClustCheat ( hitHandle,         evt, PFParticleLabelCheat );
  art::FindManyP< recob::Shower  > PFtoShowerCheat ( PFPartHandleCheat, evt, PFParticleLabelCheat );
  
  // ==================================================================
  // Run the energy reconstructions
  // ==================================================================
  
//  int nPrimaries = 0;
//  for ( size_t pf{0} ; pf < PFPartList.size() ; pf++ ) {
//    auto particle = PFPartList[pf];
//    if ( particle->IsPrimary() ) nPrimaries++;
//    const std::vector< size_t > & daughters = particle->Daughters();
//    std::cout << "\nThis is PFParticle " << pf << ".\n"
//              << "It has PDG Code " << particle->PdgCode() << ".\n"
//              << "Is this PFParticle primary? " << (particle->IsPrimary() ? "Yes" : "No") << ".\n"
//              << "Number of daughters is " << particle->NumDaughters() << ".\n"
//              << "These daughters are:\n";
//    for ( size_t i{0} ; i < daughters.size() ; i++ ) {
//      std::cout << daughters[i] << "\n";
//    }
//    std::cout << "Index of this particle is " << particle->Self() << ".\n"
//              << "The Parent of this particle is " << particle->Parent() << ".\n";
//  }
//  if ( nPrimaries > 1 ) std::cout << "\n\nEvent " << evt.id().event() << " has more than one primary. It has " << nPrimaries << " primaries\n\n";
  
  // Find a discrepancy with the MC IDE's
  art::ServiceHandle<cheat::BackTrackerService>       bt_serv;
  
  std::vector<double> IDEsums = { 0, 0, 0 };
  
  for ( size_t h{0} ; h < hitList.size() ; h++ ) {
    art::Ptr<recob::Hit> hit = hitList[h];
    const std::vector< const sim::IDE * > IDEs = bt_serv->HitToSimIDEs_Ps(hit);
    for ( size_t ide{0} ; ide < IDEs.size() ; ide++ ) {
    
      if ( hit->WireID().Plane == 0 ) IDEsums[0] += IDEs[ide]->energy/1000;
      if ( hit->WireID().Plane == 1 ) IDEsums[1] += IDEs[ide]->energy/1000;
      if ( hit->WireID().Plane == 2 ) IDEsums[2] += IDEs[ide]->energy/1000;
    }
  }
  

  MCIDE = IDEsums;
  
//  std::cout << "\n\nUnique IDE sum =     " << hitIDESum << "\n"
//            << "Non unique IDE sum = " << nonUniqueIDEsum << "\n"
//            << "Plane 0 sum =        " << plane0IDEsum << "\n"
//            << "Plane 1 sum =        " << plane1IDEsum << "\n"
//            << "Plane 2 sum =        " << plane2IDEsum << "\n";
  
  
  // Go through spacepoints and see if any near edge
  for ( size_t sp{0} ; sp < spacePointList.size() ; sp++) {
    auto spacepoint = spacePointList[sp];
    std::vector<double> xyz = { 0, 0, 0 };
    xyz[0] = spacepoint->XYZ()[0];
    xyz[1] = spacepoint->XYZ()[1];
    xyz[2] = spacepoint->XYZ()[2];
    spacePoints.push_back(xyz);
  }
  
  for ( size_t sp{0} ; sp < spacePointListCheat.size() ; sp++) {
    auto spacepoint = spacePointListCheat[sp];
    std::vector<double> xyz = { 0, 0, 0 };
    xyz[0] = spacepoint->XYZ()[0];
    xyz[1] = spacepoint->XYZ()[1];
    xyz[2] = spacepoint->XYZ()[2];
    spacePointsCheat.push_back(xyz);
  }
  
  
  // emshower
  if ( emshowerHandle.isValid() || emshowerList.size() > 0 ) {
   TrackShowerEnergyReco( emshowerToHit, emshowerList, t0, 1, 0, 0 );
  }
  else {
    std::cerr << "\n\n#####################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "emshowerHandle is not valid, or empty\n"
              << "emshowerList.size(): " << emshowerList.size() << "\n"
              << "#####################################\n\n";
  }
  
  // pandoraShower
  if ( pandoraHandle.isValid() || pandoraList.size() > 0 ) {
    TrackShowerEnergyReco( pandoraToHit, pandoraList, t0, 0, 0, 0 );
  }
  else {
    std::cerr << "\n\n####################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "pandoraHandle is not valid, or empty\n"
              << "pandoraList.size(): " << pandoraList.size() << "\n"
              << "####################################\n\n";
  }
  
  // pandoraShowerCheat
  if ( pandoraHandleCheat.isValid() || pandoraList.size() > 0 ) {
    TrackShowerEnergyReco( pandoraToHitCheat, pandoraListCheat, t0, 0, 0, 1 );
  }
  else {
    std::cerr << "\n\n####################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "pandoraHandleCheat is not valid, or empty\n"
              << "pandoraListCheat.size(): " << pandoraListCheat.size() << "\n"
              << "####################################\n\n";
  }
  
  // pandoraTrack
  if ( pandoraTrackHandle.isValid() || pandoraTrackList.size() > 0 ) {
    TrackShowerEnergyReco( pandoraTrackToHit, pandoraTrackList, t0, 0, 1, 0 );
  }
  else {
    std::cerr << "\n\n####################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "pandoraTrackHandle is not valid, or empty\n"
              << "pandoraTrackList.size(): " << pandoraTrackList.size() << "\n"
              << "####################################\n\n";
  }
  
  // PFParticles
  if ( PFPartHandle.isValid() || PFPartList.size() > 0 ) {
    PFenergyReco( PFtoClust, HitToClust, hitList, PFPartList, t0, 0 );
    PFenergyRecoByFlow ( PFtoClust, HitToClust, hitList, PFPartList, t0, 0);
  }
  else {
    std::cerr << "\n\n###################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "PFPartHandle is not valid, or empty\n"
              << "PFPartList.size(): " << PFPartList.size() << "\n"
              << "###################################\n\n";
  }
  
  // Cheated PFParticles
  if ( PFPartHandleCheat.isValid() || PFPartListCheat.size() > 0 ) {
    PFenergyReco( PFtoClustCheat, HitToClustCheat, hitList, PFPartListCheat, t0, 1 );
    PFenergyRecoByFlow ( PFtoClustCheat, HitToClustCheat, hitList, PFPartListCheat, t0, 1);
  }
  else {
    std::cerr << "\n\n########################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "PFPartHandleCheat is not valid, or empty\n"
              << "PFPartListCheat.size(): " << PFPartListCheat.size() << "\n"
              << "########################################\n\n";
  }
  
  // ==================================================================
  // Fill the Root Tree
  // ==================================================================
  
  if ( PFlowPDGPrimaryCheat.size() == 0 ) {
    std::cerr << "\n\n########################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "PFlowPDGPrimaryCheat is empty\n"
              << "########################################\n\n";
  }
  if ( PFlowPDGPrimary.size() == 0 ) {
    std::cerr << "\n\n########################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "PFlowPDGPrimary is empty\n"
              << "########################################\n\n";
  }
  
  showerTree->Fill();
  
}

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::PFenergyRecoByFlow ( art::FindManyP<recob::Cluster>             const & PFtoClust,
                                                     art::FindManyP<recob::Cluster>             const & HitToClust,
                                                     std::vector< art::Ptr<recob::Hit> >        const & hitList,
                                                     std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                                     double const & t0,
                                                     bool   const & cheat )
{

  // Associate the hits to clusters
  std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > clustersToHits;

  for ( size_t hitIndex{0} ; hitIndex < hitList.size() ; hitIndex++ ) { // hitIndex Loop

    if ( !HitToClust.at(hitIndex).empty() ) {
      const art::Ptr<recob::Cluster> cluster = HitToClust.at(hitIndex)[0];
      const art::Ptr<recob::Hit> hit = hitList[hitIndex];
      clustersToHits[cluster].push_back(hit);
    }

  } // hitIndex Loop

  // Want to follow the flow of daughter particles. So must find primary particles. Then find daughter particles, and then daughters of daughters and so on...

  // Primaries
  for ( size_t PFIndex{0} ; PFIndex < PFPartList.size() ; PFIndex++ ) {
    auto particle = PFPartList[PFIndex];
    
    if ( particle->IsPrimary() ) {
      std::vector<double> energyTotal = { 0, 0, 0 };
      auto daughters = particle->Daughters();
      
      if ( !daughters.empty() ) { // Daughters
        std::vector<double> tempEnergy = DaughterEnergy(PFPartList, PFtoClust, clustersToHits, daughters, t0);
        energyTotal[0] += tempEnergy[0];
        energyTotal[1] += tempEnergy[1];
        energyTotal[2] += tempEnergy[2];
      } // Daughters
      
      std::vector<double> totalPlaneCharge = { 0, 0, 0 };

      std::vector< art::Ptr<recob::Cluster> > allClusters = PFtoClust.at( PFIndex );

      for ( size_t clustIndex{0} ; clustIndex < allClusters.size() ; clustIndex++ ) { // clusters

        const art::Ptr<recob::Cluster> cluster = allClusters.at(clustIndex);
        auto clustIt = clustersToHits.find(cluster);

        std::vector< art::Ptr<recob::Hit> > allHits = clustIt->second;

        for ( size_t h{0} ; h < allHits.size() ; h++ ) { // h loop
          art::Ptr<recob::Hit> hit = allHits[h];
          if ( hit->WireID().Plane == 0 ) totalPlaneCharge[0] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
          if ( hit->WireID().Plane == 1 ) totalPlaneCharge[1] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
          if ( hit->WireID().Plane == 2 ) totalPlaneCharge[2] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
        } // h loop
    
      } // clusters
      
      for ( size_t i{0} ; i < totalPlaneCharge.size() ; i++ ) {
        energyTotal[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
      }
    
      if ( cheat ) {
        EfromPFlowCheat.push_back( energyTotal );
        PFlowPDGPrimaryCheat.push_back( particle->PdgCode() );
      }
      else {
        EfromPFlow.push_back( energyTotal );
        PFlowPDGPrimary.push_back( particle->PdgCode() );
      }
    
    } // Primary loop
  
  } // PFIndex loop
  
}

//---------------------------------------------------------------------------------------------------------

std::vector<double> protoDUNE::ShowerStudyAna::DaughterEnergy ( std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                                                art::FindManyP<recob::Cluster>             const & PFtoClust,
                                                                std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > const & clustersToHits,
                                                                std::vector<size_t> const & daughters,
                                                                double const & t0)
{
  std::vector<double> energyTotal = { 0, 0, 0 };
  
//  for ( std::vector<size_t>::iterator it = daughters.begin() ; it != daughters.end() ; ++it ) {
  for ( size_t i{0} ; i < daughters.size() ; i++ ) {
//    auto particle = PFPartList[*it];
    auto particle = PFPartList[ daughters[i] ];
    auto newDaughters = particle->Daughters();
    
    if ( !newDaughters.empty() ) {
      std::vector<double> tempEnergy = DaughterEnergy( PFPartList, PFtoClust, clustersToHits, newDaughters, t0);
      energyTotal[0] += tempEnergy[0];
      energyTotal[1] += tempEnergy[1];
      energyTotal[2] += tempEnergy[2];
    }
    
    std::vector<double> totalPlaneCharge = { 0, 0, 0 };

    std::vector< art::Ptr<recob::Cluster> > allClusters = PFtoClust.at( daughters[i] );

    for ( size_t clustIndex{0} ; clustIndex < allClusters.size() ; clustIndex++ ) {

      const art::Ptr<recob::Cluster> cluster = allClusters.at(clustIndex);
      auto clustIt = clustersToHits.find(cluster);

      std::vector< art::Ptr<recob::Hit> > allHits = clustIt->second;

      for ( size_t h{0} ; h < allHits.size() ; h++ ) {
        art::Ptr<recob::Hit> hit = allHits[h];
        if ( hit->WireID().Plane == 0 ) totalPlaneCharge[0] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
        if ( hit->WireID().Plane == 1 ) totalPlaneCharge[1] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
        if ( hit->WireID().Plane == 2 ) totalPlaneCharge[2] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      } // h loop
  
    }
    
    for ( size_t i{0} ; i < totalPlaneCharge.size() ; i++ ) {
      energyTotal[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
    }
    
  } // iterator loop
  
  return energyTotal;
  
}

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::PFenergyReco ( art::FindManyP<recob::Cluster>             const & PFtoClust,
                                               art::FindManyP<recob::Cluster>             const & HitToClust,
                                               std::vector< art::Ptr<recob::Hit> >        const & hitList,
                                               std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                               double const & t0,
                                               bool   const & cheat )
{

  // Associate the hits to clusters
  std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > clustersToHits;
  
  for ( size_t hitIndex{0} ; hitIndex < hitList.size() ; hitIndex++ ) { // hitIndex Loop
    
    if ( !HitToClust.at(hitIndex).empty() ) {
      const art::Ptr<recob::Cluster> cluster = HitToClust.at(hitIndex)[0];
      const art::Ptr<recob::Hit> hit = hitList[hitIndex];
      clustersToHits[cluster].push_back(hit);
    }
    
  } // hitIndex Loop
  
  std::vector<double> totE_showers = { 0, 0, 0 };
  std::vector<double> totE_tracks  = { 0, 0, 0 };
  std::vector<double> totE         = { 0, 0, 0 };
  
  // Go through the PFParticles and their clusters
  for ( size_t PFIndex{0}; PFIndex < PFPartList.size() ; PFIndex++ ) {
  
    std::vector<double> energyFromCharge = { 0, 0, 0 };
    std::vector<double> totalPlaneCharge = { 0, 0, 0 };
    
    std::vector< art::Ptr<recob::Cluster> > allClusters = PFtoClust.at( PFIndex );
    
    for ( size_t clustIndex{0} ; clustIndex < allClusters.size() ; clustIndex++ ) {
      
      const art::Ptr<recob::Cluster> cluster = allClusters.at(clustIndex);
      auto clustIt = clustersToHits.find(cluster);
      
      std::vector< art::Ptr<recob::Hit> > allHits = clustIt->second;
      
      for ( size_t h{0} ; h < allHits.size() ; h++ ) {
        art::Ptr<recob::Hit> hit = allHits[h];
        if ( hit->WireID().Plane == 0 ) totalPlaneCharge[0] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
        if ( hit->WireID().Plane == 1 ) totalPlaneCharge[1] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
        if ( hit->WireID().Plane == 2 ) totalPlaneCharge[2] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      } // h loop
      
    } // clustIndex loop
    
    for ( size_t i{0} ; i < totalPlaneCharge.size() ; i++ ) {
      if ( PFPartList[PFIndex]->PdgCode() == 11 ) totE_showers[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
      if ( PFPartList[PFIndex]->PdgCode() == 13 || PFPartList[PFIndex]->PdgCode() == 211 ) totE_tracks[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
      energyFromCharge[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
      totE[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
    }
    
    if ( cheat ) {
      if ( PFPartList[PFIndex]->PdgCode() == 11 ) PFShowersCheat.push_back( energyFromCharge );
      if ( PFPartList[PFIndex]->PdgCode() == 13 || PFPartList[PFIndex]->PdgCode() == 211 ) PFTracksCheat.push_back( energyFromCharge );
    }
    else {
      if ( PFPartList[PFIndex]->PdgCode() == 11 ) PFShowers.push_back( energyFromCharge );
      if ( PFPartList[PFIndex]->PdgCode() == 13 || PFPartList[PFIndex]->PdgCode() == 211 ) PFTracks.push_back( energyFromCharge );
    }
    
  } // PFIndex loop
  
  if ( cheat ) {
    PFShowersTotCheat = totE_showers;
    PFTracksTotCheat  = totE_tracks;
    PFTotCheat        = totE;
  }
  else {
    PFShowersTot = totE_showers;
    PFTracksTot  = totE_tracks;
    PFTot        = totE;
  }
  
}

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::TrackShowerEnergyReco ( art::FindManyP<recob::Hit>             const & TSToHit,
                                                        std::vector< art::Ptr<recob::Shower> > const & TSList,
                                                        double const & t0,
                                                        bool   const & emshower,
                                                        bool   const & track,
                                                        bool   const & cheat )
{
  
  std::vector<double> totE = { 0, 0, 0 };
  
  for ( size_t TSIndex{0} ; TSIndex < TSList.size() ; TSIndex++ ) {
  
    std::vector<double> energyFromCharge = { 0, 0, 0 };
    std::vector<double> totalPlaneCharge = { 0, 0, 0 };
  
    // Get list of hits in track/shower
    std::vector< art::Ptr<recob::Hit> > allHits = TSToHit.at( TSIndex );
    
    for ( size_t h{0} ; h < allHits.size() ; h++ ) {
      art::Ptr<recob::Hit> hit = allHits[h];
      
      if ( hit->WireID().Plane == 0 ) totalPlaneCharge[0] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      if ( hit->WireID().Plane == 1 ) totalPlaneCharge[1] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      if ( hit->WireID().Plane == 2 ) totalPlaneCharge[2] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      
    } // h looop
    
    for ( size_t i{0} ; i < totalPlaneCharge.size() ; i++ ) {
      totE[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
      energyFromCharge[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
    }
    
    if ( emshower ) {
      emShower.push_back( energyFromCharge );
    }
    else if ( track ) {
      pandoraTrack.push_back( energyFromCharge );
    }
    else {
      if ( cheat ) {
        pandoraShowerCheat.push_back( energyFromCharge );
      }
      else {
        pandoraShower.push_back( energyFromCharge );
      }
    }
  
  } // TSIndex loop
  
  if ( emshower ) {
      emShowerTotE = totE;
    }
    else if ( track ) {
      pandoraTrackTotE = totE;
    }
    else {
      if ( cheat ) {
        pandoraShowerCheatTotE = totE;
      }
      else {
        pandoraShowerTotE = totE;
      }
    }
  

}

//---------------------------------------------------------------------------------------------------------

void protoDUNE::ShowerStudyAna::TrackShowerEnergyReco ( art::FindManyP<recob::Hit>             const & TSToHit,
                                                        std::vector< art::Ptr<recob::Track> >  const & TSList,
                                                        double const & t0,
                                                        bool   const & emshower,
                                                        bool   const & track,
                                                        bool   const & cheat )
{
  
  std::vector<double> totE = { 0, 0, 0 };
  
  for ( size_t TSIndex{0} ; TSIndex < TSList.size() ; TSIndex++ ) {
  
    std::vector<double> energyFromCharge = { 0, 0, 0 };
    std::vector<double> totalPlaneCharge = { 0, 0, 0 };
  
    // Get list of hits in track/shower
    std::vector< art::Ptr<recob::Hit> > allHits = TSToHit.at( TSIndex );
    
    for ( size_t h{0} ; h < allHits.size() ; h++ ) {
      art::Ptr<recob::Hit> hit = allHits[h];
      
      if ( hit->WireID().Plane == 0 ) totalPlaneCharge[0] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      if ( hit->WireID().Plane == 1 ) totalPlaneCharge[1] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      if ( hit->WireID().Plane == 2 ) totalPlaneCharge[2] += hit->Integral() * fCaloAlg.LifetimeCorrection(hit->PeakTime(), t0);
      
    } // h looop
    
    for ( size_t i{0} ; i < totalPlaneCharge.size() ; i++ ) {
      totE[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
      energyFromCharge[i] += ( totalPlaneCharge[i] * recombCorrection * ADCtoGeV );
    }
    
    if ( emshower ) {
      emShower.push_back( energyFromCharge );
    }
    else if ( track ) {
      pandoraTrack.push_back( energyFromCharge );
    }
    else {
      if ( cheat ) {
        pandoraShowerCheat.push_back( energyFromCharge );
      }
      else {
        pandoraShower.push_back( energyFromCharge );
      }
    }
  
  } // TSIndex loop
  
  if ( emshower ) {
    emShowerTotE = totE;
  }
  else {
    if ( track ) {
      pandoraTrackTotE = totE;
    }
    else {
      if ( cheat ) {
        pandoraShowerCheatTotE = totE;
      }
      else {
        pandoraShowerTotE = totE;
      }
    }
  }
  
}

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
DEFINE_ART_MODULE(protoDUNE::ShowerStudyAna)























































