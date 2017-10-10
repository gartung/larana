// \brief: Implementation file for likelihood-based PID algorithm
// \author: Andrew Olivier aoliv23@lsu.edu

#include "LikelihoodAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

pid::LikelihoodAlg::LikelihoodAlg(fhicl::ParameterSet const& pset): fPDFFile(nullptr)
{
  //Get needed parameters from a fcl parameter set
  std::string pdfFileName = pset.get<std::string>("PDFFileName"); //name of the file that contains TH2D PDFs to be used
  std::vector<int> PDGs = pset.get<std::vector<int>>("PDGs"); //These names will be used to get PDFs from the file with pdfFileName,
                                                                                        //and each name will have a likelihood result mapped to it.  
  size_t nPlanes = pset.get<size_t>("nPlanes"); //How many planes to prepare likelihoods for
  fPDFMaps = std::vector<std::map<int,TH2D*> >(nPlanes);

  //Set up the PDF file
  std::string pdfFileFullPath;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(pdfFileName,pdfFileFullPath);
  if(pdfFileFullPath.empty())
    throw cet::exception("Likelihood File Not Found") <<"Couldn't find '" << pdfFileName << "' in FW_SEARCH_PATH";
  fPDFFile = TFile::Open(pdfFileFullPath.c_str());
  if(fPDFFile) //maybe throw if this condition is not true
  { 
    //Set up the PDF maps
    for (size_t iPlane = 0; iPlane < nPlanes; iPlane++)
    {
      for(const auto& PDG: PDGs)
      {
        std::string histName = "pdg";
        histName += std::to_string(PDG);
        histName += "_plane";
        histName += std::to_string(iPlane);
        TH2D* part_tmp = (TH2D*)(fPDFFile->Get(histName.c_str())); 
        if(!part_tmp)
        {
          throw cet::exception("LikelihoodAlg:NoPDF") << "Failed to get PDF for PDG: " << PDG << " and plane: "<<iPlane<<" from file " << fPDFFile->GetName() << ". The histogram should be in the root directory of the file and named \""<< histName << "\"\n";
        }

        part_tmp = PrepHisto(part_tmp); // turn histogram into likelihood

        if(part_tmp)
        {
          fPDFMaps[iPlane].insert(std::make_pair(PDG, part_tmp)); //make a new entry in the map of PDFs to PDG IDs
          mf::LogWarning("LikelihoodAlg") << "Added (" << PDG << ", " << part_tmp->GetName() << ") to fPDFMap for plane " << iPlane << ".\n";
        }
        else 
        {
            throw cet::exception("LikelihoodAlg:FailLoadPDF") << "Failed to turn PDF into likelihood for " << PDG << " and iPlane " << iPlane << ".\n";
        }
      } //end loop over particle names
    } // end loop over iPlane
  } //end if PDF file exists
  else // couldn't get file
  {
    throw cet::exception("LikelihoodAlg:FileError") << "Failed to open file " << pdfFileFullPath << "\n";
  }
} //end constructor

pid::LikelihoodAlg::~LikelihoodAlg()
{ 
  //I think we can and should leave this empty since all dynamic objects are owned by ROOT
} //end destructor

TH2D* pid::LikelihoodAlg::PrepHisto(const TH2D* histo)
{
  TH2D* retVal = (TH2D*)(histo->Clone());
  for(int xbin = 0; xbin < retVal->GetNbinsX(); ++xbin)
  {
    for(int ybin = 0; ybin < retVal->GetNbinsY(); ++ybin)
    {
      int bin = retVal->GetBin(xbin, ybin);
      int content = retVal->GetBinContent(bin);
      if(content == 0) retVal->SetBinContent(bin, 1e-1);
    }
  }
   double integral = retVal->Integral(0, retVal->GetNbinsX(), 0, retVal->GetNbinsY());
  if(integral != 0)
  {
    retVal->Scale(1/integral);
  }
  mf::LogWarning("LikelihoodAlg") << "In PrepHisto, got histogram named " << retVal->GetName() << " ready for likelihood.\n";
  return retVal;
}

std::map<int, double> pid::LikelihoodAlg::CalcLikelihood(const anab::Calorimetry &calo)
{
  std::map<int, double> result;
  size_t iPlane = calo.PlaneID().Plane;
  if (iPlane > 100)
  {
          mf::LogWarning("LikelihoodAlg") << "plane " << iPlane << " seems like too many planes, returning empty result.\n";
          return result;
  }
  if (iPlane >= fPDFMaps.size())
  {
    throw cet::exception("LikelihoodAlg::CalcLikelihood:InvalidPlane") << "This calo object is for plane "<< iPlane << "while we only have "<<fPDFMaps.size() <<"planes\n";
  }

  for(const auto& pdfPair: fPDFMaps[iPlane]) //calculate a likelihood value for calo for each PDF
  {
    result.insert(std::make_pair(pdfPair.first, CalcLikelihood(pdfPair, calo))); //Use private overload to calculate a likelihood value for the Calorimetry
                                                                                 //object calo using the PDF mapped to the PDG ID in pdfPair
  } //end loop over PDG-PDF map
  return result;
} //end public overload of CalcLikelihood

double pid::LikelihoodAlg::CalcLikelihood(std::pair<int, TH2D*> pdfPair, const anab::Calorimetry &calo) //calculate log likelihood for an 
                                                                                                                  //individual particle
{
  const auto resR = calo.fResidualRange;
  const auto dedx = calo.fdEdx;
  double likelihood = 0.; 
  bool goodLikelihoodPoints = false;
  TH2D* hist = pdfPair.second; //temporary pointer to PDF to be used

  TAxis* resRAxis = hist->GetXaxis();
  int resRNbins = resRAxis->GetNbins();
  double resRMin = resRAxis->GetBinLowEdge(1);
  double resRMax = resRAxis->GetBinLowEdge(resRNbins);
  TAxis* dedxAxis = hist->GetYaxis();
  int dedxNbins = dedxAxis->GetNbins();
  double dedxMin = dedxAxis->GetBinLowEdge(1);
  double dedxMax = dedxAxis->GetBinLowEdge(dedxNbins);

  //assuming that calorimetry object vectors are simultaneous
  for(size_t it = 0; it < resR.size() && it < dedx.size(); ++it)
  {
    double dedxVal = dedx[it];
    double resRVal = resR[it];
    if(resRVal < resRMin || resRVal > resRMax) continue;
    if(dedxVal < dedxMin || dedxVal > dedxMax) continue;
    double content = hist->GetBinContent(hist->FindBin(resRVal, dedxVal));
    likelihood += log(content); //log of product is sum of logs
    goodLikelihoodPoints = true;
  }
  if (!goodLikelihoodPoints)
  {
    likelihood = -999.;
  }
  return likelihood;
} //end private overload of CalcLikelihood

std::vector<int> pid::LikelihoodAlg::getPDGs() const
{
  if (fPDFMaps.size() == 0)
  {
    throw cet::exception("LikelihoodAlg::getPDGs:NoPlanes") << "No planes so can't find PDGs\n";
  }
  std::vector<int> result;
  for (const auto& pair: fPDFMaps[0])
  {
    result.push_back(pair.first);
  }
  return result;
}
