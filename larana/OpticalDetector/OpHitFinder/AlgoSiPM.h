//=============================================================================
// AlgoSiPM.h
// This is a hit finding algorithm that makes an optical hit out of
// everything above threshold, but uses only the first peak to assign hit time
//
// Gleb Sinev, Duke, 2015
// Based on AlgoTreshold.h
//=============================================================================

#ifndef ALGOSIPM_H
#define ALGOSIPM_H

#include "fhiclcpp/ParameterSet.h"

#include "PMTPulseRecoBase.h"

#include <vector>

namespace pmtana {

  class AlgoSiPM : public PMTPulseRecoBase {
    
  public:
    
    // Default constructor
    AlgoSiPM(fhicl::ParameterSet const& pset, 
             std::string         const  name = "AlgoSiPM");
    
    // Default destructor
    ~AlgoSiPM();
    
    // Implementation of PMTPulseRecoBase::Reset() method
    void Reset();
    
    // A method to set user-defined ADC threshold value
    //      void SetADCThreshold(double v) {_adc_thres = v;};
    
    // A method to set a multiplication factor to the pedestal standard 
    // deviation which is used as one of two input values to define a threshold
    //      void SetNSigma(double v) {_nsigma = v;};
    
  protected:
    
    bool RecoPulse(pmtana::Waveform_t      const&,
                   pmtana::PedestalMean_t  const&,
                   pmtana::PedestalSigma_t const&);
    
    // Use pedestal RMS from the pedestal algortihm 
    // instead of ADC counts to set thresholds
    bool _use_ped_rms;
    
    // A variable holder for a user-defined absolute threshold value
    double _adc_thres; // In ADC counts
    double _rms_thres; // In multiples of RMS
    
    // Start recording hit information after this threshold is reached
    double _2nd_adc_thres; // In ADC counts
    double _2nd_rms_thres; // In multiples of RMS
    
    // Minimum width for a hit to be recorded
    int _min_width;
    
    // Use the pedestal value from the pedestal algortihm
    bool _use_ped_algo;
    
    // Use this pedestal value if _use_ped_algo is false
    double _pedestal;

  };

}

#endif
