/*!
 * Title:   Algorithim class for dE/dx profile based PID 
 * Author:  Robert Sulej (robert.sulej@cern.ch), Dorota Stefan (dorota.stefan@cern.ch)
 *
 * Description: Calculates probabilities of particle ID from a sequence of (dE/dx; range) points,
 * used for e/gamma separation, and stopping particle detection/classification.
*/

#include "larana/ParticleIdentification/ProfilePatternPIDAlg.h"

#include "sys/stat.h"

void pid::ProfilePatternPIDAlg::reconfigure(fhicl::ParameterSet const & p)
{
	deleteNN();
	fPatternPdgs.clear();

	fPatternFile = p.get<std::string>("PatternFile");

	struct stat buff;
	if (stat(fPatternFile.c_str(), &buff) == 0) // file exists
	{
		fNNet = new pid::NNReader(fPatternFile.c_str());

		switch (fNNet->GetOutputLength())
		{
			case 2:
				fPatternPdgs.push_back(11);
				fPatternPdgs.push_back(22);
				break;
			case 4:
				fPatternPdgs.push_back(2212);
				fPatternPdgs.push_back(321);
				fPatternPdgs.push_back(211);
				fPatternPdgs.push_back(13);
				break;
			default: throw "as of today only 2 outputs (e/gamma) or 4 outputs (p/K/pi/mu)";
		}
	}
	else std::cout << "pattern file not found" << std::endl;
}

std::map< int, double > pid::ProfilePatternPIDAlg::run(
	std::vector<double> const & dedx, std::vector<double> const & range)
{
	std::map< int, double > result;
	if (!fNNet)
	{
		TVectorT<float> input(2), outsum(fNNet->GetOutputLength());
		outsum = 0.0F;
		size_t npoints = 0;
		for (size_t i = 0; i < dedx.size(); ++i)
		{
			input[0] = range[i]; input[1] = dedx[i];
			if ((input[0] > 0.0F) && (input[1] > 0.0F) && ((input[1] < 1.0e4F)))
			{
				fNNet->Run(input);
				TVectorT<float> output = fNNet->GetAllOutputs();

				bool ok = true;
				for (int j = 0; j < output.GetNoElements(); ++j)
				{
					output[j] += std::log(output[j]);
					if (isnan(output[j]) || isinf(output[j]))
					{
						ok = false; break;
					}
				}
				if (ok)
				{
					outsum += output; npoints++;
				}
			}
		}
		if (npoints)
		{
			for (int j = 0; j < outsum.GetNoElements(); ++j)
				outsum[j] = std::exp(outsum[j]);
			outsum *= 1.0F / outsum.Sum();

			for (size_t j = 0; j < fPatternPdgs.size(); ++j)
			{
				result[fPatternPdgs[j]] = outsum[j];
			}
		}
	}
	return result;
}

