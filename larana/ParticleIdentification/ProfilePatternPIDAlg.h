#ifndef PROFILEPATTERNPIDALG_H
#define PROFILEPATTERNPIDALG_H
/*!
 * Title:   Algorithim class for dE/dx profile based PID
 * Author:  Robert Sulej (robert.sulej@cern.ch), Dorota Stefan (dorota.stefan@cern.ch)
 *
 * Description: Calculates probabilities of particle ID from a sequence of (dE/dx; range) points,
 * used for e/gamma separation, and stopping particle detection/classification.
 * The general purpose is to translate dE/dx sequence of any length into a single probability value
 * (or vector if multiple PID), that can be used by any other next level MVA algorithm that combines
 * also other features of tracks, showers, etc.
 *
 * Started on Feb.8 2016.
*/
#include "fhiclcpp/ParameterSet.h"
#include "lardata/AnalysisBase/Calorimetry.h"
#include "larana/ParticleIdentification/profilepatterns/NNReader.h"

namespace pid
{
  class ProfilePatternPIDAlg;
}

class pid::ProfilePatternPIDAlg
{
public:
	ProfilePatternPIDAlg(fhicl::ParameterSet const & p) : fNNet(0) { this->reconfigure(p); }
	~ProfilePatternPIDAlg(void) { deleteNN(); }

	void reconfigure(fhicl::ParameterSet const & p);

	std::map< int, double > run(std::vector<double> const & dedx, std::vector<double> const & range);

private:
	void deleteNN(void) { if (fNNet) delete fNNet; fNNet = 0; }
	pid::NNReader* fNNet;

	std::vector< int > fPatternPdgs;
	std::string fPatternFile;
};

#endif
