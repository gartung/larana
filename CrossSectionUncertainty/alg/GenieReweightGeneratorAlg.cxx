
#include "GenieReweightGeneratorAlg.h"

namespace rwgt{

  GenieReweightGeneratorAlg::GenieReweightGeneratorAlg(){
  }

  GenieReweightGeneratorAlg::~GenieReweightGeneratorAlg(){

    // Delete all of the instances that have been created:
    for (auto & vec : reweightVector){
      for (auto & ptr : vec ) {
        if (ptr) delete ptr;
      }
    }

  }



/**
* 
* This function configures the reweighting machinery.  I am guess that,
* due to the overhead in configuring the genie reweighting, it's faster
* to make one reweight object for each weight needed.
* This will also make a "total" reweight object that will set all the
* switches on for all weighting parameters.
*/
  void GenieReweightGeneratorAlg::configureReWeight(const std::vector<reweight> & weights,
                         const std::vector<std::vector<float>>& reweightingSigmas){

    // Expand the vector to the correct number of physical knobs plus 1
    reweightVector.resize(weights.size()+1);

    // if (weights.size()+1 != reweightingSigmas.size()){
    //   std::cerr << "Error configuring the reweights, the number of weights must be"
    //             << " equal to the number of boundaries (range).\n";
    //   exit(-1);
    // }
 
    // don't forget the last vector with all of the weights
    reweightVector.back().resize(reweightingSigmas.front().size());
    for (auto & ptr : reweightVector.back()) ptr = new rwgt::NuReweight;

    for (unsigned int i_weight = 0; i_weight < reweightingSigmas.front().size(); ++i_weight)
    {
      // reweightVector.back().at(i_weight) 
      //     -> ReweightNCEL(reweightingSigmas[kNCEL][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightQEMA(reweightingSigmas[kQEMA][i_weight]);
      // reweightVector.back().at(i_weight) 
      //     -> ReweightQEVec(reweightingSigmas[kQEVec][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightResGanged(reweightingSigmas[kResGanged][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightCCRes(reweightingSigmas[kCCRes][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightNCRes(reweightingSigmas[kNCRes][i_weight]);
      // reweightVector.back().at(i_weight) 
      //     -> ReweightCoh(reweightingSigmas[kCoh][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightNonResRvp1pi(reweightingSigmas[kNonResRvp1pi][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightNonResRvbarp1pi(reweightingSigmas[kNonResRvbarp1pi][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightNonResRvp2pi(reweightingSigmas[kNonResRvp2pi][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightNonResRvbarp2pi(reweightingSigmas[kNonResRvbarp2pi][i_weight]);
      // reweightVector.back().at(i_weight) 
      //     -> ReweightResDecay(reweightingSigmas[kResDecay][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightNC(reweightingSigmas[kNC][i_weight]);
      // reweightVector.back().at(i_weight) 
      //     -> ReweightDIS(reweightingSigmas[kDIS][i_weight]);
      reweightVector.back().at(i_weight) 
          -> ReweightDISnucl(reweightingSigmas[kDISnucl][i_weight]);
      // reweightVector.back().at(i_weight) 
      //     -> ReweightAGKY(reweightingSigmas[kAGKY][i_weight]);
    }


    // loop over the physical knobs and expand to the correct number of weights
    for (unsigned int i_reweightingKnob = 0;
         i_reweightingKnob < reweightVector.size()-1; 
         i_reweightingKnob++) 
    {

      // resize this row to accomodate all of the weight points 
      reweightVector[i_reweightingKnob].resize(
            reweightingSigmas[i_reweightingKnob].size());

      for (unsigned int weight_point = 0; 
           weight_point < reweightingSigmas[i_reweightingKnob].size(); 
           weight_point++){

        // Figure out what is the value going in to this reweight
        // double stepSize = (reweightingSigmas[i_reweightingKnob] 
        //                - rangeLow[i_reweightingKnob])/(nWeights-1);
        // double reweightingValue = rangeLow[i_reweightingKnob] 
        //                         + weight_point*stepSize;

        reweightVector[i_reweightingKnob][weight_point] = new rwgt::NuReweight;

        switch (weights[i_reweightingKnob]){
          case kNCEL:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightNCEL(reweightingSigmas[kNCEL][weight_point]);
            break;
          case kQEMA:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightQEMA(reweightingSigmas[kQEMA][weight_point]);
            break;
          case kQEVec:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightQEVec(reweightingSigmas[kQEVec][weight_point]);
            break;
          case kResGanged:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightResGanged(reweightingSigmas[kResGanged][weight_point]);
            break;
          case kCCRes:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightCCRes(reweightingSigmas[kCCRes][weight_point]);
            break;
          case kNCRes:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNCRes(reweightingSigmas[kNCRes][weight_point]);
            break;
          case kCoh:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightCoh(reweightingSigmas[kCoh][weight_point]);
            break;
          case kNonResRvp1pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvp1pi(reweightingSigmas[kNonResRvp1pi][weight_point]);
            break;
          case kNonResRvbarp1pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvbarp1pi(reweightingSigmas[kNonResRvbarp1pi][weight_point]);
            break;
          case kNonResRvp2pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvp2pi(reweightingSigmas[kNonResRvp2pi][weight_point]);
            break;
          case kNonResRvbarp2pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvbarp2pi(reweightingSigmas[kNonResRvbarp2pi][weight_point]);
            break;
          case kResDecay:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightResDecay(reweightingSigmas[kResDecay][weight_point]);
            break;
          case kNC:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNC(reweightingSigmas[kNC][weight_point]);
            break;
          case kDIS:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightDIS(reweightingSigmas[kDIS][weight_point]);
            break;
          case kDISnucl:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightDISnucl(reweightingSigmas[kDISnucl][weight_point]);
            break;
          case kAGKY:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightAGKY(reweightingSigmas[kAGKY][weight_point]);
            break;
          case kNReWeights:
            break;
        }

      } //loop over nWeights
    } // loop over physical knobs




    std::cout << "\n\n\nsetup finished, running"
              << " configure on each weight......\n\n\n";

    // Tell all of the reweight drivers to configure themselves:
    for(auto & vec : reweightVector){
      for (auto & driver : vec){
        driver -> Configure();
      }
    }

    return;

  }

  void GenieReweightGeneratorAlg::prepareSigmas(int NWeights, 
                           unsigned int RandSeed,
                           std::vector<std::vector<float> > & reweightingSigmas)
  {

    TRandom rand;
    rand.SetSeed(RandSeed);

    
    reweightingSigmas.resize(kNReWeights);
    for (unsigned int i = 0; i < reweightingSigmas.size(); ++i)
    {
      reweightingSigmas[i].resize(NWeights);
      for (int j = 0; j < NWeights; j ++)
        reweightingSigmas[i][j] = rand.Gaus(0,1);
    }
    return;
  }

  void GenieReweightGeneratorAlg::parseWeights(const std::vector<std::string> & string_weights,
                              std::vector<reweight> & enum_weights){

    enum_weights.resize(string_weights.size());
    for( auto & s : string_weights){
      if (s == "NCEL") enum_weights.push_back(kNCEL);
      else if (s == "QEMA") enum_weights.push_back(kQEMA);
      else if (s == "QEVec") enum_weights.push_back(kQEVec);
      else if (s == "ResGanged") enum_weights.push_back(kResGanged);
      else if (s == "CCRes") enum_weights.push_back(kCCRes);
      else if (s == "NCRes") enum_weights.push_back(kNCRes);
      else if (s == "Coh") enum_weights.push_back(kCoh);
      else if (s == "NonResRvp1pi") enum_weights.push_back(kNonResRvp1pi);
      else if (s == "NonResRvbarp1pi") enum_weights.push_back(kNonResRvbarp1pi);
      else if (s == "NonResRvp2pi") enum_weights.push_back(kNonResRvp2pi);
      else if (s == "NonResRvbarp2pi") enum_weights.push_back(kNonResRvbarp2pi);
      else if (s == "ResDecay") enum_weights.push_back(kResDecay);
      else if (s == "NC") enum_weights.push_back(kNC);
      else if (s == "DIS") enum_weights.push_back(kDIS);
      else if (s == "DISnucl") enum_weights.push_back(kDISnucl);
      else if (s == "AGKY") enum_weights.push_back(kAGKY);
      else {
        /* REALLY NEED TO MAKE THIS CONFORM TO LARSOFT PRACTICES */
        std::cout << "The physical process you requested is not available to reweight." << std::endl;
        exit(-1);
      }
    }



  }


  void GenieReweightGeneratorAlg::calcWeight(art::Ptr<simb::MCTruth> mctruth,
                            art::Ptr<simb::GTruth > gtruth,
                            std::vector<std::vector<float>>& weights){
    // return reweight.CalcWeight(*mctruth,*gtruth);
    // if (weights.size() == 0) weights.resize(1);
    // weights.front().push_back( reweight -> CalcWeight(*mctruth,*gtruth) );

    // weights needs to be the size of the reweighting vector
    if (weights.size() != reweightVector.size())
      weights.resize(reweightVector.size());

    for (unsigned int i_weight = 0; i_weight < reweightVector.size(); i_weight ++){
      if (weights[i_weight].size() != reweightVector[i_weight].size()){
        weights[i_weight].resize(reweightVector[i_weight].size());
      }
      for (unsigned int i_reweightingKnob = 0;
           i_reweightingKnob < reweightVector[i_weight].size();
           i_reweightingKnob ++)
      {
        weights[i_weight][i_reweightingKnob] 
          = reweightVector[i_weight][i_reweightingKnob] 
            -> CalcWeight(*mctruth,*gtruth);
      }
    }
    
    return;
  }

} // end of namespace lar1nd













