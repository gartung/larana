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
#include "lardataobj/AnalysisBase/Calorimetry.h"

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
    std::string fCalorimetryModuleLabel;
    
    TTree* fEventTree;
    Int_t run;
    Int_t subrun;
    Int_t event;
    Int_t npfparticles;
    Int_t trk_id;
    Float_t start_x;
    Float_t end_x;
    Float_t start_y;
    Float_t end_y;
    Float_t start_z;
    Float_t end_z;
    Float_t max_trklen;
    Float_t mindis_0;
    Float_t mindis_1;
    Int_t origin;
    Int_t pdg;
    Float_t flash_x_diff;
    Float_t cathode_pierce_val;
    Float_t flash_z_diff;
    Float_t delta_ll;
    Float_t av_dedx_1;
    Float_t av_dedx_2;
    Float_t low_end;
    Float_t high_end;
    Float_t cosmic_score;
    Float_t true_time;
};

neutrino::NeutrinoPFParticleTagger::NeutrinoPFParticleTagger(fhicl::ParameterSet const & p)
 :mcsfitter(p.get< fhicl::ParameterSet >("fitter"))
// Initialize member data here.
{
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<anab::CosmicTag>>();
    produces< art::Assns<anab::CosmicTag,   recob::Track>>();
    produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

neutrino::NeutrinoPFParticleTagger::~NeutrinoPFParticleTagger()
{
    // Clean up dynamic memory and other resources here.
}

void neutrino::NeutrinoPFParticleTagger::produce(art::Event & evt)
{
    // Instatiate the output
    std::unique_ptr< std::vector< anab::CosmicTag > >                  cosmicTagTrackVector(       new std::vector<anab::CosmicTag>                  );
    std::unique_ptr< art::Assns<anab::CosmicTag,   recob::Track > >    assnOutCosmicTagTrack(      new art::Assns<anab::CosmicTag,   recob::Track   >);
    std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag > > assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);
    
    // Recover handle for PFParticles
    //reset();
    
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel( fPFParticleModuleLabel, pfParticleHandle);
    
    art::FindManyP<recob::Cluster> fmcp(pfParticleHandle, evt, fPFParticleModuleLabel);
    
    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrackModuleLabel,trackListHandle);
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    
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
    // const geo::TPCGeo &tpc = geom->TPC(0);
    double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;
    for (size_t i = 0; i<geom->NTPC(); ++i){
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      const geo::TPCGeo &tpc = geom->TPC(i);
      tpc.LocalToWorld(local,world);
      if (minx>world[0]-geom->DetHalfWidth(i))
	minx = world[0]-geom->DetHalfWidth(i);
      if (maxx<world[0]+geom->DetHalfWidth(i))
	maxx = world[0]+geom->DetHalfWidth(i);
      if (miny>world[1]-geom->DetHalfHeight(i))
	miny = world[1]-geom->DetHalfHeight(i);
      if (maxy<world[1]+geom->DetHalfHeight(i))
	maxy = world[1]+geom->DetHalfHeight(i);
      if (minz>world[2]-geom->DetLength(i)/2.)
	minz = world[2]-geom->DetLength(i)/2.;
      if (maxz<world[2]+geom->DetLength(i)/2.)
	maxz = world[2]+geom->DetLength(i)/2.;
    }
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event(); 
    
//    if (!pfParticleHandle.isValid())
//    {
//        evt.put( std::move(cosmicTagTrackVector) );
//        evt.put( std::move(assnOutCosmicTagTrack) );
//        return;
//    }
    
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
        
	reset();
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
          //std::vector<sim::IDE> ides;
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
              //if(pfPartIdx < max_pfparticles) origin=ii->first;
              Trackid=ii->first;
         }
     }
	
     const simb::MCParticle* particle=bt->TrackIDToParticle(Trackid);	    
     if(particle){
       if(pfPartIdx < max_pfparticles){ 
           pdg = particle->PdgCode();
	   origin = bt->TrackIDToMCTruth(Trackid)->Origin();
	   true_time = particle->T();
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
	int trkid=-1;
	//float mindis0 = FLT_MAX;
        //float mindis1 = FLT_MAX;
	//int track_vec_index=-1; 
	float del_ll=0;
	
	std::cout << "************ PF particle ID : " << int(pfPartIdx) << std::endl; 
	
	for (size_t itrk = 0; itrk<trackVec.size(); ++itrk){
	    art::Ptr<recob::Track> track = trackVec[itrk];
	    TVector3 pos,end;
	    pos = track->Vertex();
	    end = track->End();
	    
	    if(track->Length() > long_length){ 
	       long_length=track->Length();
	       //track_vec_index=itrk;
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
	
	//if(pfPartIdx < max_pfparticles){
	anab::CosmicTagID_t    tag_id   = anab::CosmicTagID_t::kNotTagged;
	   if(trkid!=-1){
	      start_x = trk_start_x;
	      end_x = trk_end_x;
	      start_y = trk_start_y;
	      end_y = trk_end_y;
	      start_z = trk_start_z;
	      end_z = trk_end_z;
	      max_trklen = long_length;
	      trk_id = trkid;
	      delta_ll = del_ll;
	      float mindis0 = FLT_MAX;
              float mindis1 = FLT_MAX;
	      
	      
	      if(trk_start_x < 0 || trk_start_x >260) cosmic_score=1;
	      if(trk_end_x < 0 || trk_end_x >260) cosmic_score=1;
	      
	      
	      if (std::abs(trk_start_y - miny)<mindis0) mindis0 = std::abs(trk_start_y- miny);
	      if (std::abs(trk_start_y - maxy)<mindis0) mindis0 = std::abs(trk_start_y - maxy);
	      if (std::abs(trk_start_z - minz)<mindis0) mindis0 = std::abs(trk_start_z - minz);
	      if (std::abs(trk_start_z - maxz)<mindis0) mindis0 = std::abs(trk_start_z - maxz);
	      if (std::abs(trk_end_y - miny)<mindis1) mindis1 = std::abs(trk_end_y - miny);
	      if (std::abs(trk_end_y - maxy)<mindis1) mindis1 = std::abs(trk_end_y - maxy);
	      if (std::abs(trk_end_z - minz)<mindis1) mindis1 = std::abs(trk_end_z - minz);
	      if (std::abs(trk_end_z - maxz)<mindis1) mindis1 = std::abs(trk_end_z - maxz);
	      
	      mindis_0 = mindis0;
	      mindis_1 = mindis1;
	      // std::cout << "****************** tpc.min Y : " << tpc.MinY() << "  *********** tpc.max Y : " << tpc.MaxY() << std::endl;
	      // std::cout << "****************** tpc min z : " << tpc.MinZ() << " ************ tpc max z : " << tpc.MaxZ() << std::endl;
	      
	      if((trk_start_y>trk_end_y?mindis0:mindis1)<30 && del_ll < -5) cosmic_score=1;
	      
	      if(mindis0 < 30 && mindis1 < 30) cosmic_score=1;
	      
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
		 
		 flash_x_diff= min_x_diff;
		 cathode_pierce_val = min_cathode_pierce;
		 flash_z_diff = min_z_diff;
		 if(/*(trk_start_y<trk_end_y?mindis0:mindis1)<30 && */min_x_diff<1) cosmic_score = 1;
		 if(/*(trk_start_y<trk_end_y?mindis0:mindis1)<30*/(mindis0 < 30 || mindis1 < 30) && min_cathode_pierce<1) cosmic_score = 1;
	      }
	      
	      else{ 
	         flash_x_diff=500;
		 cathode_pierce_val=500;
		 flash_z_diff =500;
	      }
	      
	      flash_index.clear();
	    
	      
	      /////////////////////// Calorimetry information for tracks ////////////////
	      
	      std::cout <<  "*************** Using calorimetry information first time ******************" << std::endl;
	      
	      std::vector<art::Ptr<anab::Calorimetry>> calos = fmcal.at(trkid);
	      //std::cout << "$$$$$$$$$$$$$$ Track vector ID : " << unsigned(track_vec_index) << "   " << track_vec_index << "   $$$$$$$$$$$$$$$$" << std::endl;
	      int k_plane_0_hits=-1;int k_plane_1_hits=-1;int k_plane_2_hits=-1;
	      for (size_t ical = 0; ical<calos.size(); ++ical){
	           if (!calos[ical]) continue;
	           if (!calos[ical]->PlaneID().isValid) continue;
	           int planenum = calos[ical]->PlaneID().Plane;
		   if (planenum<0||planenum>2) continue;
		   if (planenum == 0) k_plane_0_hits = int(calos[ical] -> dEdx().size());
	           else if (planenum == 1) k_plane_1_hits = int(calos[ical] -> dEdx().size());
		   else if (planenum == 2) k_plane_2_hits = int(calos[ical] -> dEdx().size());
	      }
	      
	      int k_best_planenum = -1;
	      if (k_plane_0_hits != -1 || k_plane_1_hits != -1 || k_plane_2_hits != -1){
	          if(k_plane_0_hits > k_plane_1_hits && k_plane_0_hits > k_plane_2_hits) k_best_planenum = 0;
                  else if(k_plane_1_hits > k_plane_0_hits && k_plane_1_hits > k_plane_2_hits) k_best_planenum = 1;
                  else if(k_plane_2_hits > k_plane_0_hits && k_plane_2_hits > k_plane_1_hits) k_best_planenum = 2;
                  else if(k_plane_2_hits==k_plane_0_hits && k_plane_2_hits > k_plane_1_hits) k_best_planenum = 2;
                  else if(k_plane_2_hits==k_plane_1_hits && k_plane_2_hits > k_plane_0_hits) k_best_planenum = 2;
                  else if(k_plane_1_hits==k_plane_0_hits && k_plane_1_hits > k_plane_2_hits) k_best_planenum = 0;
                  else if(k_plane_1_hits==k_plane_0_hits && k_plane_1_hits==k_plane_2_hits) k_best_planenum = 2;
	      }
	      
	      float tot_dedx_0=0;
	      float tot_dedx_1=0;
	      int n0=0;
	      int n1=0;
	      float y_hit_1=0;
	      float y_hit_last=0;
	      float res_1=0;
	      float res_last=0;
	      
	      for (size_t ical = 0; ical<calos.size(); ++ical){
	           if (!calos[ical]) continue;
	           if (!calos[ical]->PlaneID().isValid) continue;
	           int planenum = calos[ical]->PlaneID().Plane;
		   if (planenum<0||planenum>2) continue;
		   if (planenum == k_best_planenum){
		       const size_t NHits = calos[ical] -> dEdx().size();
		       for(size_t iHit = 0; iHit < NHits; ++iHit){
		           const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
			   if(iHit==0){ 
			       y_hit_1 = TrkPos.Y();
			       res_1 = (calos[ical]->ResidualRange())[iHit];
			   }
			   if(iHit==NHits-1){
			       y_hit_last = TrkPos.Y();
			       res_last = (calos[ical]->ResidualRange())[iHit];
			   } 
			   if((calos[ical]->ResidualRange())[iHit] < 3){
			       tot_dedx_0 += (calos[ical]->dEdx())[iHit];
			       n0++;
			   }
			   
			   if(calos[ical]->Range()-(calos[ical]->ResidualRange())[iHit] < 3){
			      tot_dedx_1 += (calos[ical]->dEdx())[iHit];
			      n1++;
			   }
		       }
		   }
	      }
	      
	      if(n0!=0) av_dedx_1 = float(tot_dedx_0)/n0;
	      if(n1!=0) av_dedx_2 = float(tot_dedx_1)/n1;
	      
	    //////////////////////// End of calorimetry information for tracks ////////
	    
	     if(n0!=0 && n1!=0){
	        if(y_hit_1 < y_hit_last){
		   if(res_1 < res_last){
		      low_end = float(tot_dedx_0)/n0;
		      high_end = float(tot_dedx_1)/n1;
		   }
		   
		   if(res_1 > res_last){
		     low_end = float(tot_dedx_1)/n1;
		     high_end = float(tot_dedx_0)/n0;
		   }
		}
		if(y_hit_1 > y_hit_last){
		   if(res_1 < res_last){
		      low_end = float(tot_dedx_1)/n1;
		      high_end = float(tot_dedx_0)/n0;
		   }
		   
		   if(res_1 > res_last){
		      low_end = float(tot_dedx_0)/n0;
		      high_end = float(tot_dedx_1)/n1;
		   }
		}
		
		if((start_y>end_y?mindis_0:mindis_1)<30 && (start_y<end_y?mindis_0:mindis_1)>30){
		   if(low_end >3 && high_end <3){
		     cosmic_score=1;
		   }
		}
	     }
	    
	    fEventTree->Fill();
            std::vector<float> endPt1;
            std::vector<float> endPt2;
            endPt1.push_back( trk_start_x );
            endPt1.push_back( trk_start_y );
            endPt1.push_back( trk_start_z );
            endPt2.push_back( trk_end_x );
            endPt2.push_back( trk_end_y );
            endPt2.push_back( trk_end_z );

            cosmicTagTrackVector->emplace_back( endPt1, endPt2, cosmic_score, tag_id);

            util::CreateAssn(*this, evt, *cosmicTagTrackVector, trackVec,  *assnOutCosmicTagTrack );

            util::CreateAssn(*this, evt, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);

	    }//if (trkid!=-1)
           else{
             std::vector<float> tempPt1, tempPt2;
             tempPt1.push_back(-999);
             tempPt1.push_back(-999);
             tempPt1.push_back(-999);
             tempPt2.push_back(-999);
             tempPt2.push_back(-999);
             tempPt2.push_back(-999);
             cosmicTagTrackVector->emplace_back( tempPt1, tempPt2, 0., anab::CosmicTagID_t::kNotTagged);
             util::CreateAssn(*this, evt, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
           }
	//}//if (pfPartIdx < max_pfparticles)
	
	//std::cout<< "Track vector size : " << trackVec.size() << std::endl;
	//std::cout<< "Vector index is : " << int(track_vec_index) << std::endl; 
	//std::cout<< "Lognest Track Length is "<<long_length << std::endl;
	//std::cout<< "Track start & end X : " << trk_start_x << "  " << trk_end_x << std::endl;
        
   }     
    
    //npfparticles = int(pfParticleHandle->size());
    //fEventTree->Fill();

    evt.put( std::move(cosmicTagTrackVector)      );
    evt.put( std::move(assnOutCosmicTagTrack)     );
    evt.put( std::move(assnOutCosmicTagPFParticle));

    return;

} // end of produce

void neutrino::NeutrinoPFParticleTagger::reset(){
     //npfparticles=0;
     //run = -99999;
     //subrun = -99999;
     //event = -99999;
     
     //for (int i = 0; i < max_pfparticles; i++){
          start_x = -9999;
	  end_x = -9999;
	  start_y = -9999;
	  end_y = -9999;
	  start_z = -9999;
	  end_z = -9999;
	  max_trklen = -9999;
	  mindis_0 = -9999;
	  mindis_1 = -9999;
	  origin = -9999;
	  trk_id = -9999;
	  pdg = -9999;
	  flash_x_diff = -9999;
	  cathode_pierce_val = -9999;
	  flash_z_diff = -9999;
	  delta_ll = -9999;
	  av_dedx_1 = -9999;
	  av_dedx_2 = -9999;
	  low_end = -9999;
	  high_end = -9999;
	  cosmic_score = 0;
	  true_time = -9999;
     //}
}

void neutrino::NeutrinoPFParticleTagger::beginJob()
{
 art::ServiceHandle<art::TFileService> tfs;
 
 fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
 fEventTree->Branch("run",&run,"run/I");
 fEventTree->Branch("subrun",&subrun,"subrun/I");
 fEventTree->Branch("event",&event,"event/I");
 //fEventTree->Branch("npfparticles",&npfparticles,"npfparticles/I");
 fEventTree->Branch("trk_id",&trk_id,"trk_id/I");
 fEventTree->Branch("start_x",&start_x,"start_x/F");
 fEventTree->Branch("end_x",&end_x,"end_x/F");
 fEventTree->Branch("start_y",&start_y,"start_y/F");
 fEventTree->Branch("end_y",&end_y,"end_y/F");
 fEventTree->Branch("start_z",&start_z,"start_z/F");
 fEventTree->Branch("end_z",&end_z,"end_z/F");
 fEventTree->Branch("mindis_0",&mindis_0,"mindis_0/F");
 fEventTree->Branch("mindis_1",&mindis_1,"mindis_1/F");
 fEventTree->Branch("max_trklen",&max_trklen,"max_trklen/F");
 fEventTree->Branch("origin",&origin,"origin/I");
 fEventTree->Branch("pdg",&pdg,"pdg/I");
 fEventTree->Branch("flash_x_diff",&flash_x_diff,"flash_x_diff/F");
 fEventTree->Branch("cathode_pierce_val",&cathode_pierce_val,"cathode_pierce_val/F");
 fEventTree->Branch("flash_z_diff",&flash_z_diff,"flash_z_diff/F");
 fEventTree->Branch("delta_ll",&delta_ll,"delta_ll/F");
 fEventTree->Branch("av_dedx_1",&av_dedx_1,"av_dedx_1/F");
 fEventTree->Branch("av_dedx_2",&av_dedx_2,"av_dedx_2/F");
 fEventTree->Branch("low_end",&low_end,"low_end/F");
 fEventTree->Branch("high_end",&high_end,"high_end/F");
 fEventTree->Branch("cosmic_score",&cosmic_score,"cosmic_score/F");
 fEventTree->Branch("true_time",&true_time,"true_time/F");
}

void neutrino::NeutrinoPFParticleTagger::reconfigure(fhicl::ParameterSet const & p)
{
    fPFParticleModuleLabel = p.get< std::string >("PFParticleModuleLabel");
    fTrackModuleLabel      = p.get< std::string >("TrackModuleLabel", "track");
    fFlashModuleLabel      = p.get< std::string >("FlashModuleLabel");
    fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel");
}

void neutrino::NeutrinoPFParticleTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(neutrino::NeutrinoPFParticleTagger)

//////////////////////////////////////////////////////////////////////////////////////////////////////
