// \brief: Algorithm for Tagging a reconstructed track as a pion, muon, proton, or kaon based on dE/dx versus residual range using a likelihood function
// \author: Andrew Olivier aoliv23@lsu.edu

#ifndef LIKELIHOODALG
#define LIKELIHOODALG

//ROOT includes
#include "TFile.h"
#include "TH2D.h"

//Framework includes
#include "fhiclcpp/ParameterSet.h"

//LArSoft includes
#include "lardataobj/AnalysisBase/Calorimetry.h"

//c++ includes
#include <map>

namespace pid 
{
  class LikelihoodAlg;
}

class pid::LikelihoodAlg
{
  public:
    LikelihoodAlg(fhicl::ParameterSet const& pset);
    ~LikelihoodAlg();
    std::map<int, double> CalcLikelihood(const anab::Calorimetry &calo); //Calculate the likelihood value for all particle types in fPDFMap.  Returns a map of particle name to 
                                                    //likelihood value.
    std::vector<int> getPDGs() const;

  private:
    double CalcLikelihood(std::pair<int, TH2D*> pdfPair, const anab::Calorimetry &calo); //calculate likelihood for an 
                                                                                                                      //individual particle

    TFile* fPDFFile; //File with the probability density functions to be used.  Name read from fcl file.
    std::vector<std::map<int, TH2D*> > fPDFMaps; //map of particle name to probability density function histogram for each wire plane

    double fRangeCut; //Ignore points with residual range less than this value.  Read from fcl file.  

    TH2D* PrepHisto(const TH2D* histo);
};

#endif
