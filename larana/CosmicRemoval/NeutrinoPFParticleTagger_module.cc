#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/geo.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larsim/MCCheater/BackTracker.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"
#include "TTree.h"
#include "TMath.h"

const int max_pfparticles = 1000;

namespace neutrino
{
    class NeutrinoPFParticleTagger;
}

class neutrino::NeutrinoPFParticleTagger : public art::EDProducer
{
public:
    explicit NeutrinoPFParticleTagger(fhicl::ParameterSet const & p);
    virtual ~NeutrinoPFParticleTagger();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;
    void reset();

private:
    std::string fPFParticleModuleLabel;
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    trkf::TrajectoryMCSFitter mcsfitter;
    
    TTree* fEventTree;
    Int_t run;
    Int_t subrun;
    Int_t event;
    Int_t npfparticles;
    Int_t trk_id[max_pfparticles];
    Float_t start_x[max_pfparticles];
    Float_t end_x[max_pfparticles];
    Float_t start_y[max_pfparticles];
    Float_t end_y[max_pfparticles];
    Float_t start_z[max_pfparticles];
    Float_t end_z[max_pfparticles];
    Float_t max_trklen[max_pfparticles];
    Float_t mindis_0[max_pfparticles];
    Float_t mindis_1[max_pfparticles];
    Int_t origin[max_pfparticles];
    Int_t pdg[max_pfparticles];
    Float_t flash_x_diff[max_pfparticles];
    Float_t cathode_pierce_val[max_pfparticles];
    Float_t flash_z_diff[max_pfparticles];
    Float_t delta_ll[max_pfparticles];
};

neutrino::NeutrinoPFParticleTagger::NeutrinoPFParticleTagger(fhicl::ParameterSet const & p)
 :mcsfitter(p.get< fhicl::ParameterSet >("fitter"))
// Initialize member data here.
{
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    //produces< std::vector<anab::CosmicTag>>();
    //produces< art::Assns<anab::CosmicTag,   recob::Track>>();
    //produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

neutrino::NeutrinoPFParticleTagger::~NeutrinoPFParticleTagger()
{
    // Clean up dynamic memory and other resources here.
}

void neutrino::NeutrinoPFParticleTagger::produce(art::Event & evt)
{
    // Instatiate the output
    //std::unique_ptr< std::vector< anab::CosmicTag > >                  cosmicTagTrackVector(       new std::vector<anab::CosmicTag>                  );
    //std::unique_ptr< art::Assns<anab::CosmicTag,   recob::Track > >    assnOutCosmicTagTrack(      new art::Assns<anab::CosmicTag,   recob::Track   >);
    //std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag > > assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);
    
    // Recover handle for PFParticles
    reset();
    
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel( fPFParticleModuleLabel, pfParticleHandle);
    
    art::FindManyP<recob::Cluster> fmcp(pfParticleHandle, evt, fPFParticleModuleLabel);
    
    art::ServiceHandle<cheat::BackTracker> bt;
    
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fFlashModuleLabel, flashListHandle)) art::fill_ptr_vector(flashlist, flashListHandle);
    
    std::vector<float> flash_x;
    geo::PlaneID pid(0, 0, 0);
    auto const* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    for(size_t i=0; i < flashlist.size(); i++){
        flash_x.push_back(det->ConvertTicksToX(flashlist[i]->Time()/det->SamplingRate()*1e3 + det->GetXTicksOffset(pid),pid));
    }
    
    auto const* geom = lar::providerFrom<geo::Geometry>();
    const geo::TPCGeo &tpc = geom->TPC(0);
    
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event(); 
    
    if (!pfParticleHandle.isValid())
    {
        //evt.put( std::move(cosmicTagTrackVector) );
        //evt.put( std::move(assnOutCosmicTagTrack) );
        return;
    }
    
    // Recover the handle for the tracks
    //art::Handle<std::vector<recob::Track> > trackHandle;
    //evt.getByLabel( fTrackModuleLabel, trackHandle);
    
    /*if (!trackHandle.isValid())
    {
        evt.put( std::move(cosmicTagTrackVector) );
        evt.put( std::move(assnOutCosmicTagTrack) );
        return;
    }
    
    // Recover handle for track <--> PFParticle associations
    art::Handle< art::Assns<recob::PFParticle, recob::Track> > pfPartToTrackHandle;
    evt.getByLabel(fTrackModuleLabel, pfPartToTrackHandle);*/
    
    // Recover the list of associated tracks
    art::FindManyP<recob::Track> pfPartToTrackAssns(pfParticleHandle, evt, fTrackModuleLabel);
    
    // and the hits
    //art::FindManyP<recob::Hit>  hitsSpill(trackHandle, evt, fTrackModuleLabel);
    
    // The outer loop is going to be over PFParticles
    for(size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++)
    {
        
	std::vector< art::Ptr<recob::Hit> > allHits;
      //Get all hits through associated clusters
      std::vector< art::Ptr<recob::Cluster>> allClusters = fmcp.at(pfPartIdx);
      art::FindManyP<recob::Hit> fmhcl(allClusters, evt, fPFParticleModuleLabel);
      for (size_t iclu = 0; iclu<allClusters.size(); ++iclu){
        std::vector< art::Ptr<recob::Hit>> hits = fmhcl.at(iclu);
        allHits.insert(allHits.end(), hits.begin(), hits.end());
      }
      
      std::map<int,double> trkide;
      for(size_t h = 0; h < allHits.size(); ++h){
          art::Ptr<recob::Hit> hit = allHits[h];
          std::vector<sim::IDE> ides;
          std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
    
    for(size_t e = 0; e < TrackIDs.size(); ++e){
      trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
    }
      }
      // Work out which IDE despoited the most charge in the hit if there was more than one.
      double maxe = -1;
      double tote = 0;
      int Trackid = 0;
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
          tote += ii->second;
          if ((ii->second)>maxe){
               maxe = ii->second;
              //if(pfPartIdx < max_pfparticles) origin[pfPartIdx]=ii->first;
              Trackid=ii->first;
         }
     }
	
     const simb::MCParticle* particle=bt->TrackIDToParticle(Trackid);	    
     if(particle){
       if(pfPartIdx < max_pfparticles){ 
           pdg[pfPartIdx] = particle->PdgCode();
	   origin[pfPartIdx] = bt->TrackIDToMCTruth(Trackid)->Origin();
       }
     }	
	
	art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfPartIdx);
        
        // Recover the track vector
        std::vector<art::Ptr<recob::Track> > trackVec = pfPartToTrackAssns.at(pfPartIdx);
        
        // Is there a track associated to this PFParticle?
        /*if (trackVec.empty())
        {
            // We need to make a null CosmicTag to store with this PFParticle to keep sequencing correct
            std::vector<float> tempPt1, tempPt2;
            tempPt1.push_back(-999);
            tempPt1.push_back(-999);
            tempPt1.push_back(-999);
            tempPt2.push_back(-999);
            tempPt2.push_back(-999);
            tempPt2.push_back(-999);
            cosmicTagTrackVector->emplace_back( tempPt1, tempPt2, 0., anab::CosmicTagID_t::kNotTagged);
            util::CreateAssn(*this, evt, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
            continue;
        }*/
	
	float long_length =-1;
	float trk_start_x=0;
	float trk_end_x=0;
	float trk_start_y=0;
	float trk_end_y=0;
	float trk_start_z=0;
	float trk_end_z=0;
	int trkid=0;
	float mindis0 = FLT_MAX;
        float mindis1 = FLT_MAX;
	int track_vec_index=-1; 
	float del_ll=0;
	
	std::cout << "************ PF particle ID : " << int(pfPartIdx) << std::endl; 
	
	for (size_t itrk = 0; itrk<trackVec.size(); ++itrk){
	    art::Ptr<recob::Track> track = trackVec[itrk];
	    TVector3 pos,end;
	    pos = track->Vertex();
	    end = track->End();
	    
	    if(track->Length() > long_length){ 
	       long_length=track->Length();
	       track_vec_index=itrk;
	       trk_start_x=pos.X();
	       trk_end_x=end.X();
	       trk_start_y=pos.Y();
	       trk_end_y=end.Y();
	       trk_start_z=pos.Z();
	       trk_end_z=end.Z();
	       trkid = track.key();
	       recob::MCSFitResult result = mcsfitter.fitMcs(*track);
	       if(trk_start_y > trk_end_y) del_ll=result.fwdLogLikelihood()-result.bwdLogLikelihood();
	       else del_ll=result.bwdLogLikelihood()-result.fwdLogLikelihood();
	    }
	}
	
	if(pfPartIdx < max_pfparticles){
	
	   if(track_vec_index!=-1){
	      start_x[pfPartIdx] = trk_start_x;
	      end_x[pfPartIdx] = trk_end_x;
	      start_y[pfPartIdx] = trk_start_y;
	      end_y[pfPartIdx] = trk_end_y;
	      start_z[pfPartIdx] = trk_start_z;
	      end_z[pfPartIdx] = trk_end_z;
	      max_trklen[pfPartIdx] = long_length;
	      trk_id[pfPartIdx] = trkid;
	      delta_ll[pfPartIdx] = del_ll;
	      
	      
	      if (std::abs(trk_start_y - tpc.MinY())<mindis0) mindis0 = std::abs(trk_start_y- tpc.MinY());
	      if (std::abs(trk_start_y - tpc.MaxY())<mindis0) mindis0 = std::abs(trk_start_y - tpc.MaxY());
	      if (std::abs(trk_start_z - tpc.MinZ())<mindis0) mindis0 = std::abs(trk_start_z - tpc.MinZ());
	      if (std::abs(trk_start_z - tpc.MaxZ())<mindis0) mindis0 = std::abs(trk_start_z - tpc.MaxZ());
	      if (std::abs(trk_end_y - tpc.MinY())<mindis1) mindis1 = std::abs(trk_end_y - tpc.MinY());
	      if (std::abs(trk_end_y - tpc.MaxY())<mindis1) mindis1 = std::abs(trk_end_y - tpc.MaxY());
	      if (std::abs(trk_end_z - tpc.MinZ())<mindis1) mindis1 = std::abs(trk_end_z - tpc.MinZ());
	      if (std::abs(trk_end_z - tpc.MaxZ())<mindis1) mindis1 = std::abs(trk_end_z - tpc.MaxZ());
	      
	      mindis_0[pfPartIdx] = mindis0;
	      mindis_1[pfPartIdx] = mindis1;
	      
	      float minz,minx,max_x;
	      if(trk_start_x > trk_end_x){ 
	         minz=trk_end_z;
		 minx=trk_end_x;
		 max_x=trk_start_x;
	      }
	      else if(trk_start_x < trk_end_x){ 
	              minz=trk_start_z;
		      minx=trk_start_x;
		      max_x=trk_end_x;
	      }
	      
	      else{ 
	          minz=trk_start_z;
		  minx=trk_start_x;
		  max_x=trk_start_x;
	      }
	      
	      
	      std::vector<int>flash_index;
	      
	      for(size_t i=0; i < flashlist.size(); i++){
	          float min_flashz = flashlist[i]->ZCenter()-flashlist[i]->ZWidth();
		  float max_flashz = flashlist[i]->ZCenter()+flashlist[i]->ZWidth();
		  if((minz < max_flashz) && (minz > min_flashz)) flash_index.push_back(i);
	      }
	   
	      if(flash_index.size()){
	         float min_x_diff=1e10;
		 float min_cathode_pierce=1e10;
		 float min_z_diff=0;
		 for(size_t i=0; i < flash_index.size(); i++){
		     if(TMath::Abs(flash_x[flash_index[i]]-minx) < min_x_diff){ 
		         min_x_diff=TMath::Abs(flash_x[flash_index[i]]-minx);
			 min_z_diff = TMath::Abs(minz-flashlist[flash_index[i]]->ZCenter());
		     }
		     
		     if(TMath::Abs(max_x-flash_x[flash_index[i]]-256) < min_cathode_pierce){
		        min_cathode_pierce = TMath::Abs(max_x-flash_x[flash_index[i]]-256);
		     }
		 }
		 
		 flash_x_diff[pfPartIdx]= min_x_diff;
		 cathode_pierce_val[pfPartIdx] = min_cathode_pierce;
		 flash_z_diff[pfPartIdx] = min_z_diff;
	      }
	      
	      else{ 
	         flash_x_diff[pfPartIdx]=500;
		 cathode_pierce_val[pfPartIdx]=500;
		 flash_z_diff[pfPartIdx] =500;
	      }
	      
	      flash_index.clear();
	    }
	}
	
	//std::cout<< "Track vector size : " << trackVec.size() << std::endl;
	//std::cout<< "Vector index is : " << int(track_vec_index) << std::endl; 
	//std::cout<< "Lognest Track Length is "<<long_length << std::endl;
	//std::cout<< "Track start & end X : " << trk_start_x << "  " << trk_end_x << std::endl;
        
   }     
    
    npfparticles = int(pfParticleHandle->size());
    fEventTree->Fill();
    return;

} // end of produce

void neutrino::NeutrinoPFParticleTagger::reset(){
     npfparticles=0;
     run = -99999;
     subrun = -99999;
     event = -99999;
     
     for (int i = 0; i < max_pfparticles; i++){
          start_x[i] = -9999;
	  end_x[i] = -9999;
	  start_y[i] = -9999;
	  end_y[i] = -9999;
	  start_z[i] = -9999;
	  end_z[i] = -9999;
	  max_trklen[i] = -9999;
	  mindis_0[i] = -9999;
	  mindis_1[i] = -9999;
	  origin[i] = -9999;
	  trk_id[i] = -9999;
	  pdg[i] = -9999;
	  flash_x_diff[i] = -9999;
	  cathode_pierce_val[i] = -9999;
	  flash_z_diff[i] = -9999;
	  delta_ll[i] = -9999;
     }
}

void neutrino::NeutrinoPFParticleTagger::beginJob()
{
 art::ServiceHandle<art::TFileService> tfs;
 
 fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
 fEventTree->Branch("run",&run,"run/I");
 fEventTree->Branch("subrun",&subrun,"subrun/I");
 fEventTree->Branch("event",&event,"event/I");
 fEventTree->Branch("npfparticles",&npfparticles,"npfparticles/I");
 fEventTree->Branch("trk_id",trk_id,"trk_id[npfparticles]/I");
 fEventTree->Branch("start_x",start_x,"start_x[npfparticles]/F");
 fEventTree->Branch("end_x",end_x,"end_x[npfparticles]/F");
 fEventTree->Branch("start_y",start_y,"start_y[npfparticles]/F");
 fEventTree->Branch("end_y",end_y,"end_y[npfparticles]/F");
 fEventTree->Branch("start_z",start_z,"start_z[npfparticles]/F");
 fEventTree->Branch("end_z",end_z,"end_z[npfparticles]/F");
 fEventTree->Branch("mindis_0",mindis_0,"mindis_0[npfparticles]/F");
 fEventTree->Branch("mindis_1",mindis_1,"mindis_1[npfparticles]/F");
 fEventTree->Branch("max_trklen",max_trklen,"max_trklen[npfparticles]/F");
 fEventTree->Branch("origin",origin,"origin[npfparticles]/I");
 fEventTree->Branch("pdg",pdg,"pdg[npfparticles]/I");
 fEventTree->Branch("flash_x_diff",flash_x_diff,"flash_x_diff[npfparticles]/F");
 fEventTree->Branch("cathode_pierce_val",cathode_pierce_val,"cathode_pierce_val[npfparticles]/F");
 fEventTree->Branch("flash_z_diff",flash_z_diff,"flash_z_diff[npfparticles]/F");
 fEventTree->Branch("delta_ll",delta_ll,"delta_ll[npfparticles]/F");
}

void neutrino::NeutrinoPFParticleTagger::reconfigure(fhicl::ParameterSet const & p)
{
    fPFParticleModuleLabel = p.get< std::string >("PFParticleModuleLabel");
    fTrackModuleLabel      = p.get< std::string >("TrackModuleLabel", "track");
    fFlashModuleLabel      = p.get< std::string >("FlashModuleLabel");
}

void neutrino::NeutrinoPFParticleTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(neutrino::NeutrinoPFParticleTagger)

//////////////////////////////////////////////////////////////////////////////////////////////////////
