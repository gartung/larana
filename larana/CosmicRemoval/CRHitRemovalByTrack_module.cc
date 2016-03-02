////////////////////////////////////////////////////////////////////////
//
// Class:       CRHitRemovalByTrack
// Module Type: producer
// File:        CRHitRemovalByTrack_module.cc
//
// This module produces RecoBase/Hit objects after removing those
// deemed to be due to CR muons.
//
// Configuration parameters:
//
// CosmicProducerLabels    - a list of cosmic ray producers which should be or'ed
// HitProducerLabel        - the producer of the recob::Hit objects
// TrackProducerLabel      - the producer of the recob::Track objects
// CosmicTagThresholds     - a vector of thresholds to apply to label as cosmic
// EndTickPadding          - # ticks to "pad" the end tick to account for possible
//                           uncertainty in drift velocity
//
// tjyang@fnal.gov
// Modified from CRHitRemoval_module.cc
//
////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBaseArt/HitCreator.h"
#include "lardata/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// Local functions.
namespace
{
  //----------------------------------------------------------------------------
  // Filter a collection of hits (set difference).
  // This function is copied from Track3DKalmanHit_module.cc
  //
  // Arguments:
  //
  // hits      - Hit collection from which hits should be removed.
  // used_hits - Hits to remove.
  //
  void FilterHits(std::vector<art::Ptr<recob::Hit>>& hits, std::vector<art::Ptr<recob::Hit>>& used_hits)
  {
    if(used_hits.size() > 0)
      {
	// Make sure both hit collections are sorted.
	std::stable_sort(hits.begin(), hits.end());
	std::stable_sort(used_hits.begin(), used_hits.end());
	
	// Do set difference operation.
	std::vector<art::Ptr<recob::Hit>>::iterator it = std::set_difference(hits.begin(), hits.end(), used_hits.begin(), used_hits.end(), hits.begin());
	
	// Truncate hit collection.
	hits.erase(it, hits.end());
      }
  }
}

class Propagator;

class CRHitRemovalByTrack : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit CRHitRemovalByTrack(fhicl::ParameterSet const & pset);
    virtual ~CRHitRemovalByTrack();

    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:
    // Methods
    void removeTaggedHits(const recob::PFParticle*                            pfParticle,
                          const art::Handle<std::vector<recob::PFParticle> >& pfParticleHandle,
                          const art::FindManyP<recob::Cluster>&               partToClusAssns,
                          const art::FindManyP<recob::Hit>&                   clusToHitAssns,
                          std::set<const recob::PFParticle*>&                 taggedParticles,
                          art::PtrVector<recob::Hit>&                         hitVec);

    // Fcl parameters.
    std::vector<std::string> fCosmicProducerLabels;    ///< List of cosmic tagger producers
    std::string              fHitProducerLabel;        ///< The full collection of hits
    std::string              fTrackProducerLabel;      ///< Track producer
    std::vector<double>      fCosmicTagThresholds;     ///< Thresholds for tagging
    
    int                      fEndTickPadding;          ///< Padding the end tick
    
    int                      fDetectorWidthTicks;      ///< Effective drift time in ticks
    int                      fMinTickDrift;            ///< Starting tick
    int                      fMaxTickDrift;            ///< Ending tick

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    int fNumCRRejects;    ///< Number of tracks produced.
};

DEFINE_ART_MODULE(CRHitRemovalByTrack)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CRHitRemovalByTrack::CRHitRemovalByTrack(fhicl::ParameterSet const & pset) :
  fNumEvent(0),
  fNumCRRejects(0)
{
    reconfigure(pset);

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);

    // Report.
    mf::LogInfo("CRHitRemovalByTrack") << "CRHitRemovalByTrack configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
CRHitRemovalByTrack::~CRHitRemovalByTrack()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CRHitRemovalByTrack::reconfigure(fhicl::ParameterSet const & pset)
{
    fCosmicProducerLabels    = pset.get<std::vector<std::string> >("CosmicProducerLabels");
    fHitProducerLabel        = pset.get<std::string>("HitProducerLabel");
    fTrackProducerLabel      = pset.get<std::string>("TrackProducerLabel");
    fCosmicTagThresholds     = pset.get<std::vector<double> >("CosmicTagThresholds");
    fEndTickPadding          = pset.get<int>("EndTickPadding", 50);
}

//----------------------------------------------------------------------------
/// Begin job method.
void CRHitRemovalByTrack::beginJob()
{
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* geo = lar::providerFrom<geo::Geometry>();
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    
    float samplingRate  = detp->SamplingRate();
    float driftVelocity = detp->DriftVelocity( detp->Efield(), detp->Temperature() ); // cm/us
    
    fDetectorWidthTicks = 2*geo->DetHalfWidth()/(driftVelocity*samplingRate/1000); 
    fMinTickDrift       = ts->TPCTDC2Tick(0.);
    fMaxTickDrift       = fMinTickDrift + fDetectorWidthTicks + fEndTickPadding;
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method. The goal is to produce a list of recob::Hit
/// objects which are a "clean" subset of all hits and which are believed to
/// be due to a neutrino interaction. It does this by considering input CosmicTag
/// objects, relating them to PFParticles/Tracks and removing the hits
/// associated to those objects which are believed to be Cosmic Rays.
///
void CRHitRemovalByTrack::produce(art::Event & evt)
{
    ++fNumEvent;
    
    // Start by looking up the original hits
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitProducerLabel, hitHandle);

    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> ChannelHitRawDigits
      (hitHandle, evt, fHitProducerLabel);
    art::FindOneP<recob::Wire> ChannelHitWires
      (hitHandle, evt, fHitProducerLabel);
    
    // If there are no hits then there should be no output
    if (!hitHandle.isValid()) return;

    std::vector< art::Ptr<recob::Hit> >  ChHits;
    art::fill_ptr_vector(ChHits, hitHandle);
    
    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this, evt,
      /* doWireAssns */ ChannelHitWires.isValid(),
      /* doRawDigitAssns */ ChannelHitRawDigits.isValid()
      );

    // Now recover thre remaining collections of objects in the event store that we need
    // Start with tracks
    art::Handle<std::vector<recob::Track> > trackHandle;
    evt.getByLabel(fTrackProducerLabel, trackHandle);
    
    // If no tracks then no point in continuing here
    if (!trackHandle.isValid())
    {
      for( size_t h = 0; h < ChHits.size(); h++ ){
      
	art::Ptr<recob::Wire> wire = ChannelHitWires.at(h);
	art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(h);
	
	// just copy it
	hcol.emplace_back(*ChHits[h], wire, rawdigits);      
      } // for      
      
      // put the hit collection and associations into the event
      hcol.put_into(evt);
      
      return;
    }

    art::FindManyP<recob::Hit> fmh(trackHandle, evt, fTrackProducerLabel);

    // Container to contain the "bad" hits...
    std::vector<art::Ptr<recob::Hit>> taggedHits;

    for(size_t idx = 0; idx != fCosmicProducerLabels.size(); idx++){

        art::Handle<std::vector<anab::CosmicTag> > cosmicHandle;
	std::vector<art::Ptr<anab::CosmicTag>> cosmiclist;
        if (evt.getByLabel(fCosmicProducerLabels[idx], cosmicHandle))
	  art::fill_ptr_vector(cosmiclist, cosmicHandle);
	
	for (size_t i = 0; i<cosmiclist.size(); ++i){
	  if (cosmiclist[i]->CosmicScore() > fCosmicTagThresholds[idx]){
	    
	    art::FindManyP<recob::Track> fmt(cosmicHandle, evt, fCosmicProducerLabels[idx]);
	    std::vector<art::Ptr<recob::Track>> tracks = fmt.at(i);
	    for (size_t j = 0; j<tracks.size(); ++j){
	      std::cout<<"track key "<<tracks[j].key()<<std::endl;
	      std::vector<art::Ptr<recob::Hit>> hits = fmh.at(tracks[j].key());
	      for (size_t k = 0; k<hits.size(); ++k){
		if (tracks[j].key()==18) std::cout<<hits[k]->WireID()<<std::endl;
		taggedHits.push_back(hits[k]);
	      }
	    }
	  }
	}
    }

    // Are there any tagged hits?
    std::cout<<taggedHits.size()<<std::endl;
    if (!taggedHits.empty()){
      // Remove the cosmic ray tagged hits
      std::cout<<ChHits.size()<<std::endl;
      FilterHits(ChHits, taggedHits);
      std::cout<<ChHits.size()<<std::endl;
        
        // Now make the new list of output hits
      for (const auto& hit : ChHits){
	// Check on out of time hits
	if (hit->PeakTimeMinusRMS() < fMinTickDrift) continue;
	if (hit->PeakTimePlusRMS()  > fMaxTickDrift) continue;
	
	art::Ptr<recob::Hit>::key_type hit_index = hit.key();
	art::Ptr<recob::Wire> wire = ChannelHitWires.at(hit_index);
	art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(hit_index);
	
	hcol.emplace_back(*hit, wire, rawdigits);
      }
    }
    // put the hit collection and associations into the event
    hcol.put_into(evt);
    
}

//----------------------------------------------------------------------------
/// End job method.
void CRHitRemovalByTrack::endJob()
{
    double aveCRPerEvent = fNumEvent > 0 ? double(fNumCRRejects) / double(fNumEvent) : 0.;
    
    mf::LogInfo("CRHitRemovalByTrack")
        << "CRHitRemovalByTrack statistics:\n"
        << "  Number of events = " << fNumEvent << "\n"
        << "  Number of Cosmic Rays found = " << fNumCRRejects
        << ", " << aveCRPerEvent << " average/event";
}

 
