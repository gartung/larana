//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file TrackPMTLHD_module.cc
// \author Ben Jones and Christie Chiu, MIT 2010
//         Heavily revised by Ben Jones and Wes Ketchum, 2013
//         and Eric Church, 2015



#include "art/Framework/Core/EDProducer.h"
#include "AnalysisBase/FlashMatch.h"


// ROOT includes.
#include <Rtypes.h>
#include <TMath.h>
#ifndef TrackPMTLHD_h
#define TrackPMTLHD_h 1



namespace trkf{
  class BezierTrack;
}

namespace recob{
  class OpFlash;
  class Track;

}


namespace opdet {
  
  bool TrackPMTLHD_tracksort(art::Ptr<recob::Track> t1, art::Ptr<recob::Track> t2);

  class TrackPMTLHD : public art::EDProducer{
  public:
    
    TrackPMTLHD(const fhicl::ParameterSet&);
    virtual ~TrackPMTLHD();
    
    void produce(art::Event&);
    void reconfigure(fhicl::ParameterSet const& p);
      
    std::vector<double>               GetMIPHypotheses(trkf::BezierTrack* BTrack, double XOffset=0);
    std::vector<std::vector<double> > ScanHypotheses(const art::Ptr<recob::Track> t);
    void                              PrintHypotheses(std::vector<std::vector<double> > TrackHypotheses);
    double                            GetNegLLHD(std::vector<double> signal, std::vector<double> hypothesis, double UpperLim=0);
    double                            GetMinChi2(std::vector<std::vector<double> > ScannedHypotheses, std::vector<double> FlashShape);

    
    void StoreFlashMatches(std::vector<art::Ptr<recob::Track> >& Tracks, std::vector<art::Ptr<recob::OpFlash> >& Flashes, std::vector<anab::FlashMatch>& Matches, art::Event& evt);

    
    void beginJob();
    
    
  private:
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    int         fBezierResolution;
    int         fPairingMode;
    double      fLengthCut;
    double      fPECut;

    float Box(const double &, const double&);
    float fModBoxA;
    float fModBoxB;
  };

  

}

#endif




////////////////////////////////////////////////////////////////////////
/// \file  TrackPMTLHD_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(TrackPMTLHD)

}//end namespace opdet
////////////////////////////////////////////////////////////////////////



//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file TrackPMTLHD.cxx
// \author Ben Jones and Christie Chiu, MIT 2010
//
//

// LArSoft includes
#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoBase/OpFlash.h"

// FMWK includes
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/LArProperties.h"
#include "TH2D.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// Debug flag; only used during code development.
const bool debug = true;

namespace opdet {

  //-------------------------------------------------

  TrackPMTLHD::TrackPMTLHD(fhicl::ParameterSet const& pset)
  {

    produces< std::vector<anab::FlashMatch> >();
    produces< art::Assns<recob::Track, anab::FlashMatch> >();
    produces< art::Assns<recob::OpFlash, anab::FlashMatch> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void TrackPMTLHD::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");   
    fFlashModuleLabel = pset.get<std::string>("FlashModuleLabel");
    fBezierResolution = pset.get<int>("BezierResolution");
    fLengthCut        = pset.get<double>("LengthCut");
    fPECut            = pset.get<double>("PECut");
    fPairingMode      = pset.get<int>("PairingMode");
  

  }


  //-------------------------------------------------

  void TrackPMTLHD::beginJob()
  {
  }



  //-------------------------------------------------

  TrackPMTLHD::~TrackPMTLHD()
  {
  }


  //-------------------------------------------------
  std::vector<std::vector<double> > TrackPMTLHD::ScanHypotheses(const art::Ptr<recob::Track>  track)
  {

    double MinX   = 0;    //temporary
    double MaxX   = 250;  //temporary
    size_t XSteps = 50;  //temporary

    art::ServiceHandle<geo::Geometry> geom;

    std::vector<std::vector<double> > SummedPhotons(XSteps);
    for(size_t i=0; i!=XSteps; ++i)
      SummedPhotons[i].resize(geom->NOpDets());
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    float totalLength(0.0);
    float OldVertex   = track->Vertex()[0];
    
    std::vector<bool> ValidTrajectory(XSteps, true);
    art::ServiceHandle<util::LArProperties> larp;

    // skip s=0th point, cuz we need segments.
    for (size_t s=1; s!=track->NumberTrajectoryPoints(); s++)
      {
	double MIPYield   = larp->ScintYield();
	// should get from OpDigiServices->QE(), which is a step better.
	double QE         = 0.01; 
	/*
	double dQdx       = 0.0;
	if (s<track->NumberdQdx(geo::kW)) 
	  dQdx = track->DQdxAtPoint(s,geo::kW);
	*/
	double PromptFrac = 0.25;
	// calculate number of photons as 1-NumElectrons
	//	double PromptScintYield = MIPYield * QE * larp->ModBoxCorrection(dQdx) * PromptFrac;
	double PromptScintYield = MIPYield * QE * 2.1 * PromptFrac;
	if (PromptScintYield < 0.0) 
	  PromptScintYield = 0.0;
	TVector3 p, pp; // dir
	TVector3 xyz, xyzp; //pos
	track->TrajectoryAtPoint(s,xyz,p);
	track->TrajectoryAtPoint(s-1,xyzp,pp);

       
	float LightAmount     = PromptScintYield * (xyzp-xyz).Mag();
	totalLength += (xyzp-xyz).Mag();
	// for each location in x Ben/we sum(s) photons from over the whole track (in the for over b)
	// at a given PMT. Which is just what we want for the likelihood function.
	
	for(size_t i=0; i!=XSteps; ++i)
	  {
	    if(ValidTrajectory[i])
	      {
		float NewVertex = MinX + float(i)/float(XSteps)*(MaxX-MinX);	
		xyz[0] += NewVertex-OldVertex;
		
		if((xyz[0] > MaxX) || (xyz[0] < MinX) ) ValidTrajectory[i]=false; 
		double xyzv[3]; 
		xyz.GetXYZ(&xyzv[0]);
		const std::vector<float>* PointVisibility = pvs->GetAllVisibilities(xyzv);
		
		for(size_t OpDet =0; OpDet!=PointVisibility->size();  OpDet++)
		  {

		    //		    if (PointVisibility->at(OpDet) < 0.0000000001 &&  xyz[2]>0.0 )
		      //		      std::cout<<" ScanHypothesis() check : PointVisibility is 0 for opdet " << OpDet << " at x,y,z: "  << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "." << " Distance, LightAmount are " << pvs->DistanceToOpDet(xyzv,OpDet) << ", " << LightAmount << "." << std::endl;
		    
		    SummedPhotons.at(i).at(OpDet) += PointVisibility->at(OpDet) * LightAmount;
		  }
	      }
	  } // end x scan
      } // end trajectory

    for(size_t i=0; i!=SummedPhotons.size(); ++i)
      if(!ValidTrajectory[i]) SummedPhotons.at(i).clear();

    for(size_t i=0; i!=SummedPhotons.size() && ValidTrajectory[i] ; ++i) 
      {
	std::cout<<" ScanHypothesis() Sum of photons for this track for x offset : "  << i << std::endl;
	std::cout<<" ScanHypothesis() TrackLength : "  << totalLength << std::endl;
	for(size_t OpDet =0; OpDet!=geom->NOpDets();  OpDet++)
	  {
	    std::cout<<" \t OpDet : "  << OpDet << " , total photons: " << SummedPhotons.at(i).at(OpDet) << std::endl;
	  }

      }

    return SummedPhotons;
  }

   



  //-------------------------------------------------
  // Vestigal method from TrackTimeAssoc_module.cc. We leave it for posterity, though it is not called.
  //
  // This method gets the minimum chi2 achievable by marginalizing over the track x position
  // The track dQdx is allowed to vary from the MIP value and upwards
  //
  double TrackPMTLHD::GetMinChi2(std::vector<std::vector<double> > ScannedHypotheses, std::vector<double> FlashShape)
  {
    double MinChi2  = 1.E6;
    if(FlashShape.size()==0) return MinChi2;
    
    // this is the x-scan
    for(size_t i=0; i!=ScannedHypotheses.size(); ++i) 
      // ScannedHypothesis.at(i) is  a vector of length N_PMTs  for a given x_i
      {
	if(ScannedHypotheses.at(i).size()>0)
	  {
	    double Chi2 = GetNegLLHD(FlashShape, ScannedHypotheses.at(i), MinChi2);
	    if(Chi2 < MinChi2)
	      {
		MinChi2 = Chi2;
	      }
	  }
      }
    // This is where we chuck all chi2 in the x-scan and keep only the best one.
    // I need to actually hold the whole vector or make the plot here of the full Chi2.
    return MinChi2;
  }

  //-------------------------------------------------


  // Get a hypothesis for the light collected for a bezier track
  std::vector<double> TrackPMTLHD::GetMIPHypotheses(trkf::BezierTrack* Btrack, double XOffset)
  {
    art::ServiceHandle<geo::Geometry> geom;
    std::vector<double> ReturnVector(geom->NOpDets(),0);
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    float TrackLength = Btrack->GetLength();

    double xyz[3];
    for (int b=0; b!=fBezierResolution; b++)
      {
	float s               = float(b) / float(fBezierResolution);
	float dQdx            = 2.1;    // Assume MIP value
	
	Btrack->GetTrackPoint(s,xyz);
	xyz[0]+=XOffset;
	const std::vector<float>* PointVisibility = pvs->GetAllVisibilities(xyz);
	float LightAmount = dQdx*TrackLength/float(fBezierResolution);
	
	for(size_t OpDet =0; OpDet!=PointVisibility->size();  OpDet++)
	  {
	    ReturnVector.at(OpDet)+= PointVisibility->at(OpDet) * LightAmount;
	  }
      }
    return ReturnVector;
  }

  //-------------------------------------------------


  void TrackPMTLHD::produce(art::Event& evt)
  {
    
    //    int EventID = evt.id().event();
    
    
    // Read in flashes from the event
    art::Handle< std::vector<recob::OpFlash> > flashh;
    evt.getByLabel(fFlashModuleLabel, flashh);
    std::vector<art::Ptr<recob::OpFlash> > Flashes;
    for(unsigned int i=0; i < flashh->size(); ++i)
      {
	art::Ptr<recob::OpFlash> flash(flashh,i);
        if(flash->TotalPE()>fPECut) Flashes.push_back(flash);
      }

    // Read in tracks from the event
    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);
    std::vector<art::Ptr<recob::Track> >  Tracks;
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
      	if(track->Length()>fLengthCut) Tracks.push_back(track);
      }

    std::sort(Tracks.begin(), Tracks.end(), TrackPMTLHD_tracksort);

        
    art::ServiceHandle<geo::Geometry> geom;
    size_t NOpDets = geom->NOpDets();
    
    std::map<int, bool> OnBeamFlashes;
    
    std::vector<std::vector<std::vector<double> > > TrackHypotheses;
    std::vector<std::vector<double> > FlashShapes;
       

    // For each track
    for (size_t i=0; i!=Tracks.size(); ++i)
      {
	// ScanHypothesis holds an array of N x M numbers for a given Track.
	// N is the number of samples in the scan over x; M is the number of PMTs (32).
	TrackHypotheses.push_back(ScanHypotheses(Tracks.at(i) ) );
      }
    

    for(size_t f=0; f!=Flashes.size(); ++f)
      {

	std::vector<double> ThisFlashShape(NOpDets,0);
	    
	//	if(Flashes.at(f)->InBeamFrame())
	//	  {
	for(size_t i=0; i!=NOpDets; ++i)
	  ThisFlashShape[i]=Flashes.at(f)->PE(i);
	if(Flashes.at(f)->OnBeamTime()) OnBeamFlashes[f]=true;
	//	  }
	FlashShapes.push_back(ThisFlashShape);

      }

    // This map stores the Chi2 for every match:
    //    Chi2Map[Track][Flash] = chi2
    std::map<int, std::map<int, double> > Chi2Map;

    // This map sorts the preferences of each track for each flash
    //    SortedPrefs[TrackID][Chi2] = {flashid1, flashid2, ...}
    std::vector<std::map<double, std::vector<int> > > SortedPrefs(Tracks.size());
    

    for(size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	for(size_t j=0; j!=FlashShapes.size(); ++j)	    
	  {
	    double Chi2 = GetMinChi2(TrackHypotheses.at(i), FlashShapes.at(j));
	    
	    Chi2Map[i][j]=Chi2;
	    SortedPrefs[i][Chi2].push_back(j);

	  }
      }


    // This will hold the list of matches
    std::vector<anab::FlashMatch> Matches;
	


    if(fPairingMode==0)
      {
	// In pairing mode 0, store all combinatoric matches and let the user
	//   deal with Chi2s later
	for(size_t i=0; i!=Tracks.size(); ++i)
	  for(size_t j=0; j!=Flashes.size(); ++j)
 	    Matches.push_back( anab::FlashMatch(Chi2Map[i][j], j, i, (Flashes.at(j)->OnBeamTime()>0)));
      }
    
    else if(fPairingMode==1)
      {
 	// In pairing mode 1, use the stable marriage algorithm to make a guess
 	//   at good 1<->1 pairings
	
	bool StillPairing =true;
 	std::vector<int> FlashesPaired(Flashes.size(),-1);
 	std::vector<int> TracksPaired(Tracks.size(),-1);
	
        // If we made a new match in the last round, don't stop
	while(StillPairing)
	  {
	    StillPairing=false;
	    for(size_t i=0; i!=Tracks.size(); ++i)
	      {
		// If this track still to be paired
		if(TracksPaired[i]<0)
		  {
		    // Find the flash with best remaining chi2 
		    bool MadeMatch  = false;
		    for(auto itPref = SortedPrefs[i].begin(); itPref!=SortedPrefs[i].end(); ++itPref)
		      {
			for(size_t iflash =0; iflash!=itPref->second.size(); ++iflash)
			  {
			    int FlashID           = itPref->second.at(iflash);
			    int FlashExistingPref = FlashesPaired[FlashID];
			    if(FlashExistingPref < 0)
			      {
				// This flash is available - make pairing
				TracksPaired[i]         = FlashID;
				
				// if the flash is on beam, claim to be 
				//  satisfied, but don't occupy flash
				// if flash is cosmic, claim it.
				if(!OnBeamFlashes[FlashID])
				  FlashesPaired[FlashID]  = i;
				
				StillPairing = true;
				MadeMatch    = true;
			      }
			    else
			      {
				// This flash taken - flash gets to vote
				if(Chi2Map[i][FlashID] < Chi2Map[FlashExistingPref][FlashID])
				  {
				    // If the flash prefers the new guy, switch
				    FlashesPaired[FlashID]          = i;
				    TracksPaired[i]                 = FlashID;
				    TracksPaired[FlashExistingPref] = -1;
				    MadeMatch    = true;
				    StillPairing = true;
				    break;
				  }
				// or else just roll on...			    
			      }
			  }
			if(MadeMatch)  break;
		      } // end loop over chi2s
		  } // end if unpaired track
	      } // end loop over tracks
	  } // end loop until no more pairing
	
		
	for(size_t i=0; i!=Tracks.size(); ++i)
	  {
	    if(TracksPaired[i]>0)
	      {
		int TrackID = i;
		int FlashID = TracksPaired[i];
		
		Matches.push_back( anab::FlashMatch(Chi2Map[TrackID][FlashID], FlashID, TrackID, Flashes.at(FlashID)->OnBeamTime() ));
	      }
	  }
      }


    StoreFlashMatches(Tracks, Flashes, Matches, evt);
    
  }


  //--------------------------------------------------
  // Calculate the llhd between a flash hypothesis and a track 
  //  Optional : stop counting if Chi2 bigger than UpperLim
  //
  double TrackPMTLHD::GetNegLLHD(std::vector<double> signal, std::vector<double> hypothesis, double UpperLim)
  {
       
    double sumNllhd=0;
    static const double worstNllhd(1.0E6);    

    double SignalIntegral = 0;
    double HypoIntegral = 0;

    // This renormalization should not be necessary. EC, 11-Jan-2014.
    for(size_t i=0; i!=signal.size(); ++i)
      {
	SignalIntegral+=signal.at(i);
	HypoIntegral+=hypothesis.at(i);
      }

    // If light yield indicates >1 MIP light yield, 
    //   Normalize the hypothesis to the signal size.
    //   We do not allow <1 MIP hypotheses.

    double NormFactor = SignalIntegral/HypoIntegral;
    if(NormFactor > 1) 
      {
	for(size_t i=0; i!=hypothesis.size(); ++i)
	  {
	    hypothesis.at(i) *= NormFactor;
	  }
      }
    

    for(size_t i=0; i!=signal.size(); ++i)
      {    

	// Here we calculate the probability that our expected number of photons mu, 
	// (summed over all track segments from a track translated to a given x) fluctuates 
	// to that which is observed, signal i. 
	//
	// It is important that the signal=0 tubes in i contribute ...
	// The thing that makes something infinitely improbable and thus the negllhd is huge,
	// is if we expect precisely 0 signal and yet we have a non-zero contribution in 
	// the flash on a given tube, i. This destroys all hope of this flash being matched
	// to this track.

	bool range(true);
	//	double mu = std::max(hypothesis.at(i),0.00001); // zero is bad.
	double mu = hypothesis.at(i); 
	double sumProb(0.0);

	// Add all probabilities of observed fluctuating to expectation mu over
	// a range of observed +/- some error.
	int obs = std::max(std::floor(signal.at(i) - sqrt(hypothesis.at(i))), 0.0);
	//int obs = 0;
	while (range)
	  {
	    double nllhd( TMath::Poisson(obs, mu) );
	    sumProb += nllhd;
	    obs++;
	    range=false;
	    if (obs < std::ceil(signal.at(i) + sqrt(hypothesis.at(i))) ) range=true;
	    //	    if (obs <= std::floor(signal.at(i)) ) range=true;
	  }

	    double nllhd( -std::log( sumProb ) );
	    nllhd = (std::isinf(nllhd) ? worstNllhd : nllhd);
	    //	    if (nllhd == worstNllhd)
	    //	      std::cout << "GetNegLLHD: i, mu, signal: " << i <<  ", " << mu <<", " << signal.at(i) << "." << std::endl;
	    sumNllhd+=nllhd;

      }  // sum on each observed signal on PMT_i

    return sumNllhd;
  }
  

  //--------------------------------------------------

  void TrackPMTLHD::PrintHypotheses(std::vector<std::vector<double> > TrackHypotheses)
  {
    // List the light hypotheses per track, per PMT
    for (size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	mf::LogVerbatim("TrackPMTLHD")<< "Visbility for track " << i <<std::endl;

	for(size_t j=0; j!=TrackHypotheses.at(i).size(); ++j)
	  {
	    mf::LogVerbatim("TrackPMTLHD") << "Signal at PMT " << j << ", "  << TrackHypotheses.at(i).at(j)<<std::endl;
	  }
      }
    
  }



  //--------------------------------------------------

  void TrackPMTLHD::StoreFlashMatches(std::vector<art::Ptr<recob::Track> > & Tracks, std::vector<art::Ptr<recob::OpFlash> > & Flashes, std::vector<anab::FlashMatch>& Matches, art::Event& evt)
  {
    std::unique_ptr< std::vector<anab::FlashMatch> > flash_matches ( new std::vector<anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::Track, anab::FlashMatch > > assn_track( new art::Assns<recob::Track, anab::FlashMatch>);
    std::unique_ptr< art::Assns<recob::OpFlash, anab::FlashMatch > > assn_flash( new art::Assns<recob::OpFlash, anab::FlashMatch>);

    for(size_t i=0; i!=Matches.size(); ++i)
      {
	flash_matches->push_back(Matches.at(i));
	
	util::CreateAssn(*this, evt, *(flash_matches.get()), Tracks.at(Matches.at(i).SubjectID()), *(assn_track.get()), i); 

	util::CreateAssn(*this, evt, *(flash_matches.get()), Flashes.at(Matches.at(i).FlashID()), *(assn_flash.get()), i); 
      }
    

    evt.put(std::move(flash_matches));
    evt.put(std::move(assn_track));
    evt.put(std::move(assn_flash));
  }


  //-----------------------------------
  // Helper function for performing track length sort
  bool TrackPMTLHD_tracksort(art::Ptr<recob::Track> t1, art::Ptr<recob::Track> t2)
  {
    return (t1->Length()>t2->Length());
    
  }

  
}


