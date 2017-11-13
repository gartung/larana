////////////////////////////////////////////////////////////////////////
// Class:       CosmicPFParticleTagger
// Module Type: producer
// File:        CosmicPFParticleTagger_module.cc
//              This module checks timing and TPC volume boundaries as a
//              way to tag potential cosmic rays
//              This particular module uses PFParticles as input and handles
//              special cases associated with them.
//              This module started life as CosmicTrackTagger_module, written
//              by Sarah Lockwitz, and small alterations made to handle the
//              PFParticle input
//
// Generated at Wed Sep 17 19:17:00 2014 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicPFParticleTagger
////////////////////////////////////////////////////////////////////////

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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"

// Including cheater to help identify if particles being tagged are cosmics
#include "larsim/MCCheater/BackTracker.h"


namespace cosmic
{
    class CosmicPFParticleTagger;
    class SpacePoint;
    class Track;
}


class cosmic::CosmicPFParticleTagger : public art::EDProducer
{
public:
    explicit CosmicPFParticleTagger(fhicl::ParameterSet const & p);
    virtual ~CosmicPFParticleTagger();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

private:
    std::string fGenieGenModuleLabel;
    std::string fPFParticleModuleLabel;
    std::string fTrackModuleLabel;
    std::vector< art::Ptr<simb::MCTruth> > fMC;
    int         fEndTickPadding;
    int         fDetectorWidthTicks;
    int         fMinTickDrift, fMaxTickDrift;
    float       fTPCXBoundary, fTPCYBoundary, fTPCZBoundary;
    float       fDetYMin, fDetYMax, fDetHeight, fDetHalfHeight;
    float       fDetXMin, fDetXMax, fDetWidth, fDetHalfWidth;
    float       fDetZMin, fDetZMax, fDetLength;

};


cosmic::CosmicPFParticleTagger::CosmicPFParticleTagger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<anab::CosmicTag>>();
    produces< art::Assns<anab::CosmicTag,   recob::Track>>();
    produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

cosmic::CosmicPFParticleTagger::~CosmicPFParticleTagger()
{
    // Clean up dynamic memory and other resources here.
}

void cosmic::CosmicPFParticleTagger::produce(art::Event & evt)
{              
    // Instatiate the output
    std::unique_ptr< std::vector< anab::CosmicTag > >                  cosmicTagTrackVector(       new std::vector<anab::CosmicTag>                  );
    std::unique_ptr< art::Assns<anab::CosmicTag,   recob::Track > >    assnOutCosmicTagTrack(      new art::Assns<anab::CosmicTag,   recob::Track   >);
    std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag > > assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);
    
    // Recover handle for PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel( fPFParticleModuleLabel, pfParticleHandle);
    
    if (!pfParticleHandle.isValid())
    {
        evt.put( std::move(cosmicTagTrackVector) );
        evt.put( std::move(assnOutCosmicTagTrack) );
        return;
    }
    
    // Recover the handle for the tracks
    art::Handle<std::vector<recob::Track> > trackHandle;
    evt.getByLabel( fTrackModuleLabel, trackHandle);
    
    if (!trackHandle.isValid())
    {
        evt.put( std::move(cosmicTagTrackVector) );
        evt.put( std::move(assnOutCosmicTagTrack) );
        return;
    }
    
    // Recover handle for track <--> PFParticle associations
    art::Handle< art::Assns<recob::PFParticle, recob::Track> > pfPartToTrackHandle;
    evt.getByLabel(fTrackModuleLabel, pfPartToTrackHandle);
    
    // Recover the list of associated tracks
    art::FindManyP<recob::Track> pfPartToTrackAssns(pfParticleHandle, evt, fTrackModuleLabel);
    
    // and the hits
    art::FindManyP<recob::Hit>  hitsSpill(trackHandle, evt, fTrackModuleLabel);
    
    // The outer loop is going to be over PFParticles
    for(size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++)
    {
        art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfPartIdx);
        
        // Recover the track vector
        std::vector<art::Ptr<recob::Track> > trackVec = pfPartToTrackAssns.at(pfPartIdx);
        
        // Is there a track associated to this PFParticle?
        if (trackVec.empty())
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
        }
        
        // Start the tagging process...
        int                    cosmicScore =  0;
        anab::CosmicTagID_t    tag_id   = anab::CosmicTagID_t::kNotTagged;
        art::Ptr<recob::Track> track1   = trackVec.front();
        
        std::vector<art::Ptr<recob::Hit> > hitVec = hitsSpill.at(track1.key());
        
        // Recover track end points
        TVector3 vertexPosition  = track1->Vertex();
        TVector3 vertexDirection = track1->VertexDirection();
        TVector3 endPosition     = track1->End();
        
        // In principle there is one track associated to a PFParticle... but with current
        // technology it can happen that a PFParticle is broken into multiple tracks. Our
        // aim here is to find the maximum extents of all the tracks which have been
        // associated to the single PFParticle
        if (trackVec.size() > 1)
        {
            for(size_t trackIdx = 1; trackIdx < trackVec.size(); trackIdx++)
            {
                art::Ptr<recob::Track> track(trackVec[trackIdx]);
                
                TVector3 trackStart = track->Vertex();
                TVector3 trackEnd   = track->End();
                
                // Arc length possibilities for start of track
                double arcLStartToStart = (trackStart - vertexPosition).Dot(vertexDirection);
                double arcLStartToEnd   = (trackEnd   - vertexPosition).Dot(vertexDirection);
                
                if (arcLStartToStart < 0. || arcLStartToEnd < 0.)
                {
                    if (arcLStartToStart < arcLStartToEnd) vertexPosition = trackStart;
                    else                                   vertexPosition = trackEnd;
                }
                
                // Arc length possibilities for end of track
                double arcLEndToStart = (trackStart - endPosition).Dot(vertexDirection);
                double arcLEndToEnd   = (trackEnd   - endPosition).Dot(vertexDirection);
                
                if (arcLEndToStart > 0. || arcLEndToEnd > 0.)
                {
                    if (arcLEndToStart > arcLEndToEnd) endPosition = trackStart;
                    else                               endPosition = trackEnd;
                }
                
                // add the hits from this track to the collection
                hitVec.insert(hitVec.end(), hitsSpill.at(track.key()).begin(), hitsSpill.at(track.key()).end());
            }
        }

        // "Track" end points in easily readable form
        float trackEndPt1_X = vertexPosition [0];
        float trackEndPt1_Y = vertexPosition [1];
        float trackEndPt1_Z = vertexPosition [2];
        float trackEndPt2_X = endPosition[0];
        float trackEndPt2_Y = endPosition[1];
        float trackEndPt2_Z = endPosition[2];
        
        /////////////////////////////////////
        // Check that all hits on particle are "in time"
        /////////////////////////////////////
        for ( unsigned int p = 0; p < hitVec.size(); p++)
        {
            int peakLessRms = hitVec[p]->PeakTimeMinusRMS();
            int peakPlusRms = hitVec[p]->PeakTimePlusRMS();
            
            //if( hitVec[p]->PeakTimeMinusRMS() < fMinTickDrift || hitVec[p]->PeakTimePlusRMS() > fMaxTickDrift)
            if( peakLessRms < fMinTickDrift || peakPlusRms > fMaxTickDrift)
            {
                cosmicScore = 1;
                tag_id   = anab::CosmicTagID_t::kOutsideDrift_Partial;
                break;     // If one hit is out of time it must be a cosmic ray
            }
        }
        
        /////////////////////////////////
        // Now check the TPC boundaries:
        /////////////////////////////////
        if(cosmicScore==0 )
        {
          /* Below is an alternate version to the current Cosmic taggin method.
           * In this alternate method, we assign a detector face to each point, so either +-x +-y or +-z,
           * and from this we then decide whether or not to remove the particle.
           * This provides us with the power to decide whether or not to remove particles that are associated with
           * the upstream Z face (as well as any other face we might consider - maybe particles with cosmicScores of 0.5 associated
           * with the top Y face>), as with proton dune we're not using a neutrino beam, so interactions are very likely to start
           * right near the upstream Z face.
           */
          /* Create vectors of flags associating track end points to a face. 
           * -1 for the negative/smaller coordinate face, +1 for the opposite
           * 0 for if point is within the detector entirely
           */
          std::vector<int> vertexFace = {0,0,0};
          std::vector<int> endFace = {0,0,0};
          
          std::vector<float> fDetMins = {fDetXMin,fDetYMin,fDetZMin};
          std::vector<float> fDetMaxs = {fDetXMax,fDetYMax,fDetZMax};
          std::vector<float> fDetBuffers = {fTPCXBoundary,fTPCYBoundary,fTPCZBoundary};
          
          // Loop through the x,y,z positions of the points, and determine if it should be associated to a detector face
          for ( size_t i{0} ; i <  3 ; i++ ) {
            if ( vertexPosition[i] < (fDetMins[i] + fDetBuffers[i]) || vertexPosition[i] > (fDetMaxs[i] - fDetBuffers[i]) ) {
              vertexFace[i] = vertexPosition[i] < (fDetMaxs[i] + fDetMins[i])/2 ? -1 : 1;
            }
            if ( endPosition[i] < (fDetMins[i] + fDetBuffers[i]) || endPosition[i] > (fDetMaxs[i] - fDetBuffers[i]) ) {
              endFace[i] = endPosition[i] < (fDetMaxs[i] + fDetMins[i])/2 ? -1 : 1;
            }
          }
          
          // Calculate the sum of the abs values of each face associations
          int vertexFaceTot = std::abs(vertexFace[0]) + std::abs(vertexFace[1]) + std::abs(vertexFace[2]);
          int endFaceTot = std::abs(endFace[0]) + std::abs(endFace[1]) + std::abs(endFace[2]);
          
          /* Depending on the size of the detBuffer values, it's possible to get a point that is associated to two faces.
           * To solve this, we want to remove one of the associations. And we do it by prioritizing the Z association over
           * the Y, and the Y over the X 
           */
          if ( vertexFaceTot > 1 ) {
            if ( std::abs(vertexFace[0]) == 1 && std::abs(vertexFace[1]) == 1 ) {
              vertexFace[0] = 0;
            }
            else if ( std::abs(vertexFace[0]) == 1 && std::abs(vertexFace[2]) == 1 ) {
              vertexFace[0] = 0;
            }
            if ( std::abs(vertexFace[1]) == 1 && std::abs(vertexFace[2]) == 1 ) {
              vertexFace[1] = 0;
            }
          }
          if ( endFaceTot > 1 ) {
            if ( std::abs(endFace[0]) == 1 && std::abs(endFace[1]) == 1 ) {
              endFace[0] = 0;
            }
            else if ( std::abs(endFace[0]) == 1 && std::abs(endFace[2]) == 1 ) {
              endFace[0] = 0;
            }
            if ( std::abs(endFace[1]) == 1 && std::abs(endFace[2]) == 1 ) {
              endFace[1] = 0;
            }
          }        
          
          // Recalculate the sum of the abs values of each face association
          vertexFaceTot = std::abs(vertexFace[0]) + std::abs(vertexFace[1]) + std::abs(vertexFace[2]);
          endFaceTot = std::abs(endFace[0]) + std::abs(endFace[1]) + std::abs(endFace[2]);
          
          // Next we need to determine what the state of the PFParticle is.
          
          // First check those that are associated to two faces
          if ( vertexFaceTot + endFaceTot == 2 ) {
            
            if ( vertexFace[0] != 0 || endFace[0] != 0 ){
              // A particle that's going from one side of the detector out the other in X
              if ( vertexFace[0] == -endFace[0] ) {
                tag_id = anab::CosmicTagID_t::kGeometry_XX;
                cosmicScore = 1;
              }
              
              // A particle that goes from an X face to a Y face (or vice versa)
              else if ( (std::abs(vertexFace[0]) == std::abs(endFace[1])) || (std::abs(vertexFace[1]) == std::abs(endFace[0])) ) {
                tag_id = anab::CosmicTagID_t::kGeometry_XY;
                cosmicScore = 1;
              }
              
              /* A particle that is either going from X to the downstream Z face (or vice versa) or going from X to the upstream
               * Z face (or vice versa)
               */
              else if ( (std::abs(vertexFace[0]) == std::abs(endFace[2])) || (std::abs(vertexFace[2]) == std::abs(endFace[0])) ) {
                tag_id = (vertexFace[2] == 1 || endFace[2] == 1) ? anab::CosmicTagID_t::kGeometry_XDZ : anab::CosmicTagID_t::kGeometry_XUZ;
                cosmicScore = (vertexFace[2] == 1 || endFace[2] == 1) ? 1 : 0.3 ;
              }
              
              // A particle that travels parallel to the associated face. But is it completely outside the detector, or in buffer zone?
              else if ( vertexFace[0] == endFace[0] ){
                tag_id = ( (vertexPosition[0] > fDetMins[0] && vertexPosition[0] < (fDetMins[0] + fDetBuffers[0])) || (vertexPosition[0] < fDetMaxs[0] && vertexPosition[0] > (fDetMaxs[0] - fDetBuffers[0])) ) ? anab::CosmicTagID_t::kGeometry_XP : anab::CosmicTagID_t::kGeometry_XXO;
                cosmicScore = ( (vertexPosition[0] > fDetMins[0] && vertexPosition[0] < (fDetMins[0] + fDetBuffers[0])) || (vertexPosition[0] < fDetMaxs[0] && vertexPosition[0] > (fDetMaxs[0] - fDetBuffers[0])) ) ? 0.4 : 1;
              }
            }
            
            else if ( vertexFace[1] != 0 || endFace[1] != 0 ) {
              // A particle that's going from one side of the detector out the other in Y
              if ( vertexFace[1] == -endFace[1] ) {
                tag_id = anab::CosmicTagID_t::kGeometry_YY;
                cosmicScore = 1;
              }
              
              /* A particle that is either going from Y to the downstream Z face (or vice versa) or going from Y to the upstream
               * Z face (or vice versa)
               */
              else if ( (std::abs(vertexFace[1]) == std::abs(endFace[2])) || (std::abs(vertexFace[2]) == std::abs(endFace[1])) ) { 
                tag_id = ( vertexFace[2] == 1 || endFace[2] == 1 ) ? anab::CosmicTagID_t::kGeometry_YDZ : anab::CosmicTagID_t::kGeometry_YUZ;
                cosmicScore = ( vertexFace[2] == 1 || endFace[2] == 1 ) ? 1 : 0.3;
              }
              
              // A particle that travels parallel to the associated face. But is it completely outside the detector, or in buffer zone?
              else if ( vertexFace[1] == endFace[1] ) {
                tag_id = ( (vertexPosition[1] > fDetMins[1] && vertexPosition[1] < (fDetMins[1] + fDetBuffers[1])) || (vertexPosition[1] < fDetMaxs[1] && vertexPosition[1] > (fDetMaxs[1] - fDetBuffers[1])) ) ? anab::CosmicTagID_t::kGeometry_YP : anab::CosmicTagID_t::kGeometry_YYO;
                cosmicScore = ( (vertexPosition[1] > fDetMins[1] && vertexPosition[1] < (fDetMins[1] + fDetBuffers[1])) || (vertexPosition[1] < fDetMaxs[1] && vertexPosition[1] > (fDetMaxs[1] - fDetBuffers[1])) ) ? 0.4 : 1;
              }
            }
            
            else if ( vertexFace[2] != 0 || endFace[2] != 0 ) {
              // A particle that's going from one side of the detector out the other in Z
              if  ( vertexFace[2] == -endFace[2] ) {
                tag_id = anab::CosmicTagID_t::kGeometry_ZZ;
                cosmicScore = 0.3;
              }
              
              // A particle that travels parallel to the associated face. But is it completely outside the detector, or in buffer zone?
              else if ( vertexFace[2] == endFace[2] ) {
                tag_id = ( (vertexPosition[2] > fDetMins[2] && vertexPosition[2] < (fDetMins[2] + fDetBuffers[2])) || (vertexPosition[2] < fDetMaxs[2] && vertexPosition[2] > (fDetMaxs[2] - fDetBuffers[2])) ) ? anab::CosmicTagID_t::kGeometry_ZP : anab::CosmicTagID_t::kGeometry_ZZO;
                cosmicScore = ( (vertexPosition[2] > fDetMins[2] && vertexPosition[2] < (fDetMins[2] + fDetBuffers[2])) || (vertexPosition[2] < fDetMaxs[2] && vertexPosition[2] > (fDetMaxs[2] - fDetBuffers[2])) ) ? 0.4 : 1;
              }
            }
          }
          
          else if ( vertexFaceTot + endFaceTot == 1 ) {
            // A particle with one end associated to one X face, but the other is witihin the detector
            if ( vertexFace[0] != 0 || endFace[0] != 0 ) {
              tag_id = anab::CosmicTagID_t::kGeometry_X;
              cosmicScore = 0.5;
            }
            // A particle with one end associated to one Y face, but the other is witihin the detector
            else if ( vertexFace[1] != 0 || endFace[1] != 0 ) {
              tag_id = anab::CosmicTagID_t::kGeometry_Y;
              cosmicScore = 0.5;
            }
            // A particle with one end associated to one Z face, but the other is witihin the detector.
            // But is the Z face the upstream or downstream face?
            else if ( vertexFace[2] != 0 || endFace[2] != 0 ) {
              tag_id = ( vertexFace[2] == 1 || endFace[2] == 1 ) ? anab::CosmicTagID_t::kGeometry_ZD : anab::CosmicTagID_t::kGeometry_ZU;
              cosmicScore = ( vertexFace[2] == 1 || endFace[2] == 1 ) ? 0.35 : 0.3;
            }
          }
          
          else {
            cosmicScore = 0;
          }
        }
        
        std::vector<float> endPt1;
        std::vector<float> endPt2;
        endPt1.push_back( trackEndPt1_X );
        endPt1.push_back( trackEndPt1_Y );
        endPt1.push_back( trackEndPt1_Z );
        endPt2.push_back( trackEndPt2_X );
        endPt2.push_back( trackEndPt2_Y );
        endPt2.push_back( trackEndPt2_Z );
        
        // Loop through the tracks resulting from this PFParticle and mark them
        cosmicTagTrackVector->emplace_back( endPt1, endPt2, cosmicScore, tag_id);
        
        util::CreateAssn(*this, evt, *cosmicTagTrackVector, trackVec, *assnOutCosmicTagTrack );
        
        // Don't forget the association to the PFParticle
        util::CreateAssn(*this, evt, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
    }
    
    evt.put( std::move(cosmicTagTrackVector)      );
    evt.put( std::move(assnOutCosmicTagTrack)     );
    evt.put( std::move(assnOutCosmicTagPFParticle));
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void cosmic::CosmicPFParticleTagger::beginJob()
{
}

void cosmic::CosmicPFParticleTagger::reconfigure(fhicl::ParameterSet const & p)
{
    // Implementation of optional member function here.
    ////////  fSptalg  = new cosmic::SpacePointAlg(p.get<fhicl::ParameterSet>("SpacePointAlg"));
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* geo = lar::providerFrom<geo::Geometry>();
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    
    size_t nTPCs { geo->TotalNTPC() }; // Get the number of TPCs
    std::vector<double> XVec, YVec, ZVec;
    /*
     * This loop goes over each TPC and gets the max and min x,y,z values for that TPC and puts them in the vectors above.
     * This is only necessary because the geometry function DetHalfWidth(), and the like, only works on a single TPC, so 
     * doesn't do what it says on the tin (i.e. it doens't return the half width of the detector).
    */
    //std::vector<int> goodTPCs = {1,2,5,6,9,10};
    for ( size_t i{0} ; i < nTPCs ; i++ ) {
    //for (int i : goodTPCs ){
      XVec.push_back( geo->TPC(i).MaxX() );
      XVec.push_back( geo->TPC(i).MinX() );
      
      YVec.push_back( geo->TPC(i).MaxY() );
      YVec.push_back( geo->TPC(i).MinY() );
      
      ZVec.push_back( geo->TPC(i).MaxZ() );
      ZVec.push_back( geo->TPC(i).MinZ() );
    }

    fDetYMin = *std::min_element( YVec.begin(), YVec.end() );
    fDetYMax = *std::max_element( YVec.begin(), YVec.end() );
    fDetHeight = fDetYMax - fDetYMin;
    fDetHalfHeight = fDetHeight/2;

    fDetXMin = *std::min_element( XVec.begin(), XVec.end() );
    fDetXMax = *std::max_element( XVec.begin(), XVec.end() );
    fDetWidth = fDetXMax - fDetXMin;
    fDetHalfWidth = fDetWidth/2;

    fDetZMin = *std::min_element( ZVec.begin(), ZVec.end() );
    fDetZMax = *std::max_element( ZVec.begin(), ZVec.end() );
    fDetLength = fDetZMax - fDetZMin;
    
    float fSamplingRate = detp->SamplingRate();

    fGenieGenModuleLabel   = p.get< std::string >("GenieGenModuleLabel");
    fPFParticleModuleLabel = p.get< std::string >("PFParticleModuleLabel");
    fTrackModuleLabel      = p.get< std::string >("TrackModuleLabel", "track");
    fEndTickPadding        = p.get<    int      >("EndTickPadding",   50);     // Fudge the TPC edge in ticks...

    fTPCXBoundary = p.get< float >("TPCXBoundary", 5);
    fTPCYBoundary = p.get< float >("TPCYBoundary", 5);
    fTPCZBoundary = p.get< float >("TPCZBoundary", 5);

    const double driftVelocity = detp->DriftVelocity( detp->Efield(), detp->Temperature() ); // cm/us

    //std::cerr << "Drift velocity is " << driftVelocity << " cm/us.  Sampling rate is: "<< fSamplingRate << " detector width: " <<  2*geo->DetHalfWidth() << std::endl;
    fDetectorWidthTicks = fDetWidth/(driftVelocity*fSamplingRate/1000);
    fMinTickDrift = ts->TPCTDC2Tick(0.);
    fMaxTickDrift = fMinTickDrift + fDetectorWidthTicks + fEndTickPadding;
}

void cosmic::CosmicPFParticleTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(cosmic::CosmicPFParticleTagger)

