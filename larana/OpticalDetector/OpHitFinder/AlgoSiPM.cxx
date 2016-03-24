//=======================
// AlgoSiPM.cxx
//=======================

#include "AlgoSiPM.h"

#ifndef ALGOSIPM_CXX
#define ALGOSIPM_CXX

namespace pmtana {

  //---------------------------------------------------------------------------
  AlgoSiPM::AlgoSiPM(fhicl::ParameterSet const& pset,std::string const name)
    : PMTPulseRecoBase(name)
  {

    _use_ped_rms   = pset.get< bool  >("UsePedRMS");
    if (_use_ped_rms) {
      _rms_thres     = pset.get< float >("RMSThreshold"      );
      _2nd_rms_thres = pset.get< float >("SecondRMSThreshold");
    }
    else {
      _adc_thres     = pset.get< float >("ADCThreshold"      );
      _2nd_adc_thres = pset.get< float >("SecondADCThreshold");
    }
    _min_width     = pset.get< float >("MinWidth"          );
    _use_ped_algo  = pset.get< bool  >("UsePedAlgo"        );

    if (!_use_ped_algo) _pedestal = pset.get< float >("Pedestal");

    Reset();

  }

  //---------------------------------------------------------------------------
  AlgoSiPM::~AlgoSiPM()
  {}

  //---------------------------------------------------------------------------
  void AlgoSiPM::Reset()
  {

    PMTPulseRecoBase::Reset();

  }

  //---------------------------------------------------------------------------
  bool AlgoSiPM::RecoPulse(pmtana::Waveform_t      const& wf,
                           pmtana::PedestalMean_t  const& ped_mean,
                           pmtana::PedestalSigma_t const& ped_rms)
  {
    
    bool   fire           = false;
    bool   first_found    = false;
    bool   record_hit     = false;
    int    counter        = 0;
    double const pedestal = _use_ped_algo  ? 
                            // Use output of the first pedestal algorithm
                            ped_mean.at(0) :
                            _pedestal;
    double threshold      = _use_ped_rms             ?
                            // Use output of the first pedestal algorithm
                            _rms_thres*ped_rms.at(0) :
                            _adc_thres;
    threshold            += pedestal;
    double pre_threshold  = _use_ped_rms             ?
                            // Use output of the first pedestal algorithm
                            _2nd_rms_thres*ped_rms.at(0) :
                            _2nd_adc_thres;
    pre_threshold        += pedestal;

    Reset();

    for (short const &value : wf) {

      if (!fire && (double(value) >= pre_threshold)) {
        
        // Found a new pulse
        fire           = true;
        first_found    = false;
        record_hit     = false;
        _pulse.t_start = counter;

      }

      if (fire && (double(value) < pre_threshold)) {

        // Found the end of a pulse
        fire = false;
        _pulse.t_end = counter - 1;
        if (record_hit && ((_pulse.t_end - _pulse.t_start) >= _min_width)) 
        {
          _pulse_v.push_back(_pulse);
          record_hit = false;
        }
        _pulse.reset_param();

      }

      if (fire) {

        // We want to record the hit only if _adc_thres is reached
        if (!record_hit && (double(value) >= threshold)) record_hit = true;

        // Add this ADC count to the integral
        _pulse.area += (double(value) - double(pedestal));

        if (!first_found && 
            (_pulse.peak < (double(value) - double(pedestal)))) {

          // Found a new maximum
          _pulse.peak  = (double(value) - double(pedestal));
          _pulse.t_max = counter;

        }
        else if (!first_found)
          // Found the first peak
          first_found = true;

      }

      counter++;
    
    }

    if (fire) {

      // Take care of a pulse that did not finish within the readout window
      fire = false;
      _pulse.t_end = counter - 1;
      if (record_hit && ((_pulse.t_end - _pulse.t_start) >= _min_width)) 
      {
        _pulse_v.push_back(_pulse);
        record_hit = false;
      }
      _pulse.reset_param();

    }

    return true;

  }

}

#endif
