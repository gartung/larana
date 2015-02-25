#ifndef GENIE_REWEIGHT_GENERATOR_ALG_H
#define GENIE_REWEIGHT_GENERATOR_ALG_H 


#include "NuReweight/art/NuReweight.h" //GENIEReweight.h"
// #include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/GTruth.h"
// #include "SimulationBase/MCFlux.h"
/// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include <memory>
#include <iostream>

namespace rwgt{

  enum reweight {kNCEL, kQEMA, kQEVec, kResGanged, kCCRes, kNCRes, 
        kCoh, kNonResRvp1pi, kNonResRvbarp1pi, kNonResRvp2pi,
        kNonResRvbarp2pi, kResDecay, kNC, kDIS, kDISnucl, 
        kAGKY, kNReWeights};
  
  class GenieReweightGeneratorAlg
  {
  public:

    GenieReweightGeneratorAlg();
    ~GenieReweightGeneratorAlg();
  
    /**
     * @brief configures the reweights
     * @details this function does the under the hood configuration and setup of
     * all the reweighting tools with the calculated sigmas
     * 
     * @param weights the input weights
     * @param reweightingSigmas the input sigmas
     */
    void configureReWeight(const std::vector<reweight> & weights,
                           const std::vector<std::vector<float>>& reweightingSigmas);


    /**
     * @brief returns the array of weights for each particular event defined by mctruth, gtruth
     * @details The function takes mctruth and gtruth then uses a utility in nutools 
     * to convert to genie format.  It calls functions in nutools to run the actual reweighting
     * but they in turn are just wrappers over genie's functionality.  the return vector, weights,
     * is a 2 by 2 vector that always has dimensions of first index = physical parameter, second 
     * index = knobs.  The first index is always longer by 1 than the number of requested weights
     * to include a total weight.
     * 
     * @param weights [description]
     */
    void calcWeight(art::Ptr<simb::MCTruth> mctruth,
                    art::Ptr<simb::GTruth > gtruth,
                    std::vector<std::vector<float>>& weights);

    /**
     * @brief calculate sigmas to use in reweighting
     * @details Calculates a deviation from 1 for each reweighting type
     * and for each of the weights.  Allows to set the RNG seed to get reproducability.
     * 
     * @param NWeights number of weights
     * @param RandSeed the seed for the random number generator
     * @param reweightingSigmas return of 2D array of deviations (drawn from 1 sigma gaussians)
     */
    void prepareSigmas(int NWeights, unsigned int RandSeed,
                       std::vector<std::vector<float> > & reweightingSigmas);


    /**
     * @brief compares requested reweights (from .fcl) to available reweights
     * @details Loops over the entries in string_weights and does string comparison
     * to available weights.  For matches, places that weight's enum in the enum_weights
     * vector.  The vec<vec<float>> for output will store weights in an order determined
     * by this function, which is in turn determined by the order in the fcl file.
     * 
     * @param string_weights string of weights from fcl file, input
     * @param enum_weights vector of enum reweight, output, return by reference.
     */
    void parseWeights(const std::vector<std::string> & string_weights, 
                            std::vector<reweight> & enum_weights);

  private:


    // The reweighting utility class:
    std::vector<std::vector<rwgt::NuReweight *> > reweightVector;

  };
}

#endif
