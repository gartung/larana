////////////////////////////////////////////////////////////////////////
// Class:       ProfilePID
// Module Type: analyzer
// File:        ProfilePID_module.cc
//
// Author:      Robert Sulej, Feb.8 2016
//
// Examples of using PID algorithm based on dE/dx patterns only. The
// general purpose is to translate dE/dx sequence of any length into a
// single probability value (or vector if multiple PID), that can be used
// by any other next level MVA algorithm that combines also other features
// of tracks, showers, etc.
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "larsim/MCCheater/BackTracker.h"
#include "larsim/Simulation/ParticleList.h"
#include "SimulationBase/MCParticle.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/TrackHitMeta.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/PFParticle.h"

#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "larana/ParticleIdentification/ProfilePatternPIDAlg.h"

#include <functional>
#include <fstream>

#include "TTree.h"

namespace pid
{
	struct bHitInfo;
	struct bIndexLess;
	class ProfilePID;
}

struct pid::bHitInfo
{
	bHitInfo(size_t i, double x, double e) :
		Index(i), dE(e), dx(x)
	{ }
	size_t Index;
	double dE, dx;
};

struct pid::bIndexLess :
	public std::binary_function<const bHitInfo &, const bHitInfo &, bool>
{
public:
	bool operator() (const bHitInfo & h1, const bHitInfo & h2) { return h1.Index < h2.Index; }
};

class pid::ProfilePID : public art::EDAnalyzer {
public:
	explicit ProfilePID(fhicl::ParameterSet const & p);

	ProfilePID(ProfilePID const &) = delete;
	ProfilePID(ProfilePID &&) = delete;
	ProfilePID & operator = (ProfilePID const &) = delete;
	ProfilePID & operator = (ProfilePID &&) = delete;

	void analyze(art::Event const & evt) override;

	void reconfigure(fhicl::ParameterSet const & p) override;
	void beginJob() override;
	void endJob() override;

private:
	void writeTrackInfo(size_t tidx, const std::vector< double > & dEdx, const std::vector< double > & range);
	bool prepareEvent(art::Event const & evt);

	void make_dEdx(std::vector< double > & dEdx, std::vector< double > & range,
		const std::vector< pid::bHitInfo > & hits, double dvtx, double rmax) const;

	simb::MCParticle const * getMCParticle(const std::vector< art::Ptr<recob::Hit> >& trkHits, 
					double& clean, bool& stopping) const;

	art::Handle< std::vector<recob::Track> > fTrkListHandle;
	std::map< size_t, std::vector< pid::bHitInfo >[3] > fTrk2InfoMap; // hits info sorted by views
	std::map< size_t, double > fTrk2VtxDistMap; // first hit to start/end vertex distances

	int fRunNum, fEventNum, fBestView, fTrkIdx;
	double fdEdx, fRange;

	std::ofstream fPatternOutFile; // dQ/dx tree saved in Mode 1
	TTree* fPatternsTree;          // dQ/dx file saved in Mode 0

//  Parameters:
	std::string           fTrackModuleLabel;
	calo::CalorimetryAlg  fCalorimetryAlg;
	ProfilePatternPIDAlg  fProfilePIDAlg;

	std::string           fPatternOutFileName;

	int                   fDirection;
	double                fMinDx;
	double                fMaxRange;

	int                   fPdg;
	int                   fMode;
};


pid::ProfilePID::ProfilePID(fhicl::ParameterSet const & p) : EDAnalyzer(p),
	fCalorimetryAlg(p.get< fhicl::ParameterSet >("CalorimetryAlg")),
	fProfilePIDAlg(p.get<fhicl::ParameterSet>("ProfilePIDAlg"))
{
	this->reconfigure(p);
}

void pid::ProfilePID::reconfigure(fhicl::ParameterSet const & p)
{
	fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
	fCalorimetryAlg.reconfigure(p.get< fhicl::ParameterSet >("CalorimetryAlg"));
	fProfilePIDAlg.reconfigure(p.get<fhicl::ParameterSet>("ProfilePIDAlg"));

	fPatternOutFileName = p.get<std::string>("PatternOutFile");

	fDirection = p.get<int>("Direction");
	fMinDx = p.get<double>("MinDx");
	fMaxRange = p.get<double>("MaxRange");

	fPdg = p.get<int>("Pdg");
	fMode = p.get<int>("Mode");
}

void pid::ProfilePID::beginJob()
{
	if (fMode == 1) fPatternOutFile.open(fPatternOutFileName, std::ofstream::out);

	if (fMode == 0)
	{
		art::ServiceHandle<art::TFileService> tfs;
		fPatternsTree = tfs->make<TTree>("ProfilePatterns", "dE/dx info");
		fPatternsTree->Branch("fRunNum", &fRunNum, "fRunNum/I");
		fPatternsTree->Branch("fEventNum", &fEventNum, "fEventNum/I");
		fPatternsTree->Branch("fBestView", &fBestView, "fBestView/I");
		fPatternsTree->Branch("fTrkIdx", &fTrkIdx, "fTrkIdx/I");
		fPatternsTree->Branch("fdEdx", &fdEdx, "fdEdx/D");
		fPatternsTree->Branch("fRange", &fRange, "fRange/D");
	}
}

void pid::ProfilePID::endJob()
{
	if (fMode == 1) fPatternOutFile.close();
}

void pid::ProfilePID::writeTrackInfo(size_t tidx, const std::vector< double > & dEdx, const std::vector< double > & range)
{
	switch (fMode)
	{
		case 0: // save dE/dx to TTree
			for (size_t i = 0; i < dEdx.size(); ++i)
			{
				fTrkIdx = tidx; fdEdx = dEdx[i]; fRange = range[i];
				fPatternsTree->Fill();
			}
			break;

		case 1: // save dE/dx to ascii file
			for (size_t i = 0; i < dEdx.size(); ++i)
			{
				fPatternOutFile
					<< fRunNum << " " << fEventNum << " " << fBestView << " "
					<< tidx << " " << dEdx[i] << " " << range[i]
					<< std::endl;
			}
			break;

		default: return; // save data for training pattern preparation only in modes 0 and 1
	}
}

bool pid::ProfilePID::prepareEvent(art::Event const & evt)
{
	fRunNum = evt.run(); fEventNum = evt.event();
	fBestView = -1; fTrkIdx = -1;
	fdEdx = -1.0; fRange = -1.0;

	art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
	bool hasPfp = evt.getByLabel(fTrackModuleLabel, pfpListHandle);

	fTrk2InfoMap.clear();
	fTrk2VtxDistMap.clear();
	if (evt.getByLabel(fTrackModuleLabel, fTrkListHandle))
	{
		art::FindManyP< recob::Hit, recob::TrackHitMeta > hitFromTrk(fTrkListHandle, evt, fTrackModuleLabel);
		art::FindManyP< recob::PFParticle > pfpFromTrk(fTrkListHandle, evt, fTrackModuleLabel);
		art::FindManyP< recob::Vertex > vtxFromPfp(pfpListHandle, evt, fTrackModuleLabel);
		if (hitFromTrk.size())
		{
			for (size_t t = 0; t < hitFromTrk.size(); t++)
			{
				double dvtx = 0.0;
				if (hasPfp)
				{
					if (fDirection) // try to get the starting vertex (the way of searching may change...)
					{
						auto const & trk = (*fTrkListHandle)[t];
						auto pfps = pfpFromTrk.at(t);
						if (!pfps.empty())
						{
							auto vtxs = vtxFromPfp.at(pfps.front().key());
							if (!vtxs.empty())
							{
								double vtxpos[3];
								vtxs.front()->XYZ(vtxpos);
								TVector3 vtx3d(vtxpos[0], vtxpos[1], vtxpos[2]);
								dvtx = (vtx3d - trk.Vertex()).Mag();

								std::cout << "pfp vtx:" << vtx3d.X() << " " << vtx3d.Y() << " " << vtx3d.Z() << std::endl;
								std::cout << "trk vtx:" << trk.Vertex().X() << " " << trk.Vertex().Y() << " " << trk.Vertex().Z() << std::endl;
								std::cout << "dvtx:" << dvtx << std::endl;
							}
						}
					}
					else
					{
						//std::cout << "fDir:" << fDirection << " hasPfp:" << hasPfp << std::endl;
					}
				}

				auto vhit = hitFromTrk.at(t);
				auto vmeta = hitFromTrk.data(t);
				if ((fMode == 0) || (fMode == 1)) // select tracks by matched MC truth
				{
					bool stopping;
					double clean;
					auto const & mcParticle = *getMCParticle(vhit, clean, stopping);
					if (!stopping || (mcParticle.PdgCode() != fPdg) || (clean < 0.85)) continue;
				}

				fTrk2VtxDistMap[t] = dvtx;
				for (size_t h = 0; h < vhit.size(); ++h)
				{
					size_t view = vhit[h]->WireID().Plane;
					size_t idx = vmeta[h]->Index();
					double tdrift = vhit[h]->PeakTime();
					double dx = vmeta[h]->Dx();

					double dq = fCalorimetryAlg.ElectronsFromADCArea(vhit[h]->Integral(), view);
					dq *= fCalorimetryAlg.LifetimeCorrection(tdrift); // *** note: T0 = 0

					fTrk2InfoMap[t][view].emplace_back(idx, dx, dq);
				}
				for (auto & hits : fTrk2InfoMap[t]) std::sort(hits.begin(), hits.end(), pid::bIndexLess());
			}
		}
	}
	if (fTrk2InfoMap.empty()) return false;
	else return true;
}

void pid::ProfilePID::analyze(art::Event const & evt)
{
	if (!prepareEvent(evt)) return;

	std::cout << "ProfilePID: selected " << fTrk2InfoMap.size() << " tracks" << std::endl;
	for (auto const & trkEntry : fTrk2InfoMap)
	{
		size_t t = trkEntry.first;
		auto const & info = trkEntry.second;
		//auto const & trk = (*fTrkListHandle)[t];

		size_t nhits, maxhits = 0;
		for (size_t v = 0; v < 3; ++v)
		{
			nhits = info[v].size();
			if (nhits > maxhits) { maxhits = nhits; fBestView = v; }
		}
		if (fBestView < 0) continue;

		std::vector< double > dEdx, range;
		make_dEdx(dEdx, range, info[fBestView], fTrk2VtxDistMap[t], fMaxRange);

		if (fMode == 2) // apply PID
		{
			// map of PDG - probability(PID)
			std::map< int, double > pidOutput = fProfilePIDAlg.run(dEdx, range);
			std::cout << "track: " << t << " view:" << fBestView << std::endl;
			for (auto const & pid : pidOutput)
				std::cout << "   pid:" << pid.first << " p=" << pid.second << std::endl;
		}
		else writeTrackInfo(t, dEdx, range); // save training patterns
	}
}

void pid::ProfilePID::make_dEdx(std::vector< double > & dEdx, std::vector< double > & range,
	const std::vector< pid::bHitInfo > & hits, double dvtx, double rmax) const
{
	dEdx.clear(); range.clear();

	int i0, i1, di;
	if (fDirection) { i0 = 0; i1 = hits.size(); di = 1; }
	else { i0 = hits.size() - 1; i1 = -1; di = -1; }

	double de, dx, r0, r1 = 0.0, r = dvtx;
	while ((i0 != i1) && (r < rmax))
	{
		dx = 0.0; de = 0.0;
		while ((i0 != i1) && (dx <= fMinDx))
		{
			de += hits[i0].dE;
			dx += hits[i0].dx;
			i0 += di;
		}

		r0 = r1;
		r1 += dx;
		r = 0.5 * (r0 + r1);

		if ((de > 0.0) && (dx > 0.0) && (r < rmax))
		{
			dEdx.push_back(de/dx);
			range.push_back(r);
		}
	}
}

simb::MCParticle const * pid::ProfilePID::getMCParticle(
	const std::vector< art::Ptr<recob::Hit> >& trkHits, double& clean, bool& stopping) const
{
	int id = 0;
	art::ServiceHandle< cheat::BackTracker > bt;

	std::map<int, double> trkide;
	for (auto const & hit : trkHits)
		for(auto const & ide : bt->HitToTrackID(hit))
		{
			trkide[ide.trackID] += ide.energy;
		}
	
	double maxe = -1., tote = 0.;
	for (auto const & entry : trkide)
	{
		tote += entry.second;
		if ((entry.second) > maxe)
		{
			maxe = entry.second;
			id = entry.first;
		}
	}

	if (tote > 0.) clean = maxe / tote;
	else clean = 0.;

	simb::MCParticle const * particle = bt->TrackIDToParticle(id);
	double ek = particle->EndE() - particle->Mass();

	if ((ek < 0.001) || (particle->EndProcess() == "FastScintillation")) stopping = true;
	else stopping = false;

	std::cout << particle->PdgCode() << " ek:" << ek
		<< " process: " << particle->EndProcess()
		<< " clean:" << clean << " stopping:" << stopping
		<< std::endl;

	return particle;
}

DEFINE_ART_MODULE(pid::ProfilePID)
