#ifndef OPFLASHALG_H
#define OPFLASHALG_H
/*!
 * Title:   OpFlash Algorithims
 * Author:  Ben Jones, MIT (Edited by wketchum@lanl.gov)
 *
 * Description:
 * These are the algorithms used by OpFlashFinder to produce flashes.
 */

#include <functional>
#include "Simulation/BeamGateInfo.h"
#include "OpticalDetectorData/OpticalTypes.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "OpticalDetector/AlgoThreshold.h"
#include "OpticalDetector/PulseRecoManager.h"
#include "RecoBase/OpHit.h"
#include "RecoBase/OpFlash.h"

namespace opdet{

  void GetTriggerTime(std::vector<const sim::BeamGateInfo*> const&,
		      double const&,
		      double const&,
		      optdata::TimeSlice_t const&,
		      unsigned int&, unsigned short&);

  void RunFlashFinder(std::vector<optdata::FIFOChannel> const&,
		      std::vector<recob::OpHit>&,
		      std::vector<recob::OpFlash>&,
		      std::vector< std::vector<int> >&,
		      optdata::TimeSlice_t const&,
		      int const&,
		      pmtana::PulseRecoManager const&,
		      pmtana::AlgoThreshold const&,
		      std::map<int,int> const&,
		      unsigned int const&,
		      unsigned int const&,
		      float const&,
		      float const&,
		      float const&,
		      unsigned int const&,
		      unsigned int const&,
		      std::vector<double> const&);
  
  void ProcessFrame(unsigned short,
		    std::vector<const optdata::FIFOChannel*> const&,
		    std::vector<recob::OpHit>&,
		    std::vector<recob::OpFlash>&,
		    optdata::TimeSlice_t const&,
		    int const&,
		    pmtana::PulseRecoManager const&,
		    pmtana::AlgoThreshold const&,
		    std::map<int,int> const&,
		    unsigned int const&,
		    float const&,
		    float const&,
		    float const&,
		    unsigned int const&,
		    unsigned int const&,
		    std::vector<double> const&);

  void ConstructHits(int const&,
		     uint32_t const&,
		     unsigned short const&,
		     pmtana::AlgoThreshold const&,
		     std::vector<recob::OpHit>&,
		     optdata::TimeSlice_t const&,
		     int const&,
		     float const&,
		     float const&,
		     unsigned int const&,
		     unsigned int const&,
		     double const&,
		     std::vector<double> &,
		     std::vector<double> &,
		     std::vector< std::vector<int> > &,
		     std::vector< std::vector<int> > &,
		     std::vector<int> &,
		     std::vector<int> &);


  void AssignHitsToFlash( std::vector<int> const&,
			  std::vector<int> const&,
			  std::vector<double> const&,
			  std::vector<double> const&,
			  std::vector< std::vector<int> > const&,
			  std::vector< std::vector<int> > const&,
			  size_t const&,
			  std::vector<recob::OpHit> const&,
			  std::vector< std::vector<int> >&,
			  float const&);

  void RefineHitsToFlash(std::vector< std::vector<int> > const&,
			 std::vector<recob::OpHit> const&,
			 std::vector< std::vector<int> >&,
			 float const&,
			 float const&);


}//end opdet namespace

#endif
