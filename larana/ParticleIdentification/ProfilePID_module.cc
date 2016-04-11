////////////////////////////////////////////////////////////////////////
// Class:       ProfilePID
// Module Type: analyzer
// File:        ProfilePID_module.cc
//
// Author:      Robert Sulej, Feb.8 2016 and dorota.stefan@cern.ch
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
#include "larsim/Simulation/LArG4Parameters.h"
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
#include <cmath>

#include "TTree.h"

namespace pid
{
	struct bHitInfo;
	struct bIndexLess;
	class ProfilePID;
}

struct pid::bHitInfo
{
	bHitInfo(size_t i, double x, double e, int w) :
		Index(i), dE(e), dx(x), wire(w)
	{ }
	size_t Index;
	double dE, dx;
	int wire;
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
	void beginRun(const art::Run& run) override;
	void endJob() override;

private:
	void writeTrackInfo(size_t tidx, const std::vector< double > & dEdx, const std::vector< double > & range);
	bool prepareEvent(art::Event const & evt);

	bool has(const std::vector<int>& v, int i) const;

	bool evselect(art::Event const & evt) const;

	void make_dEdx(std::vector< double > & dEdx, std::vector< double > & range,
		const std::vector< pid::bHitInfo > & hits, double dvtx, double rmax) const;

	//simb::MCParticle const * getMCParticle(const std::vector< art::Ptr<recob::Hit> >& trkHits, 
	//				double& clean, bool& stopping) const;

	art::Handle< std::vector<recob::Track> > fTrkListHandle;
	std::map< size_t, std::vector< pid::bHitInfo >[3] > fTrk2InfoMap; // hits info sorted by views
	std::map< size_t, double > fTrk2VtxDistMap; // first hit to start/end vertex distances

	int fRunNum, fEventNum, fBestView, fTrkIdx;
	double fdEdx, fRange;

	double fElectronsToGeV;

	std::ofstream fPatternOutFile; // dQ/dx tree saved in Mode 1
	TTree* fPatternsTree;          // dQ/dx file saved in Mode 0

//  Parameters:
	std::string           fTrackModuleLabel;
	std::string						fSimulationLabel;
	calo::CalorimetryAlg  fCalorimetryAlg;
	ProfilePatternPIDAlg  fProfilePIDAlg;

	std::string           fPatternOutFileName;

	int                   fDirection;
	int										fFlip;
	double                fMinDx;
	double                fMaxRange;

	std::vector<int>      fPdg;
	int                   fMode;
	bool 									fStopping;
	bool									fDecaying;
};


pid::ProfilePID::ProfilePID(fhicl::ParameterSet const & p) : EDAnalyzer(p),
	fCalorimetryAlg(p.get< fhicl::ParameterSet >("CalorimetryAlg")),
	fProfilePIDAlg(p.get<fhicl::ParameterSet>("ProfilePIDAlg"))
{
	this->reconfigure(p);
}

void pid::ProfilePID::reconfigure(fhicl::ParameterSet const & p)
{
	fSimulationLabel = p.get<std::string>("SimulationLabel");
	fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
	fCalorimetryAlg.reconfigure(p.get< fhicl::ParameterSet >("CalorimetryAlg"));
	fProfilePIDAlg.reconfigure(p.get<fhicl::ParameterSet>("ProfilePIDAlg"));

	fPatternOutFileName = p.get<std::string>("PatternOutFile");

	fDirection = p.get<int>("Direction");
	fMinDx = p.get<double>("MinDx");
	fMaxRange = p.get<double>("MaxRange");

	fPdg = p.get< std::vector<int> >("Pdg");
	fMode = p.get<int>("Mode");
	fStopping = p.get<bool>("Stopping");
	fDecaying = p.get<bool>("Decaying");

	fFlip = 0;
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

void pid::ProfilePID::beginRun(const art::Run&)
{
	art::ServiceHandle<sim::LArG4Parameters> larParameters;
	fElectronsToGeV = 1./larParameters->GeVToElectrons();
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

bool pid::ProfilePID::evselect(art::Event const & evt) const
{
	// select particle based on monte carlo information
	std::vector< art::Ptr<simb::MCParticle> > simlist;

	art::Handle< std::vector<simb::MCParticle> > mcparticleHandle;
	if (evt.getByLabel(fSimulationLabel, mcparticleHandle))
		art::fill_ptr_vector(simlist, mcparticleHandle);

	std::map< int, const simb::MCParticle* > particleMap;
	for (auto const& particle : simlist)
	{
		particleMap[particle->TrackId()] = &*particle;
	}

	int simpdg = particleMap.begin()->second->PdgCode();
	if (fStopping)
	{
		if ((simpdg == 2212) 
			&& has(fPdg, simpdg) 
			&& (particleMap.size() == 1) 
			&& (particleMap.begin()->second->EndProcess() != "ProtonInelastic")) 	
		{
			return true; 
		}
	}

	if (fDecaying)
	{
		bool kaonmodedcy = false;
		if ((simpdg == 321) 
			&& has(fPdg, simpdg)
			&& (particleMap.begin()->second->EndProcess() == "Decay")
			&& (particleMap.begin()->second->NumberDaughters() == 2)
			&& (particleMap.size() == 6))
		{ 	
			size_t count = 0;
			for (auto const & p : particleMap)
			{
				if ((p.second->PdgCode() == 321) ||
						(p.second->PdgCode() == -13) ||
						(p.second->PdgCode() == 14))
				{
					kaonmodedcy = true;
				}

				count++;
				if (count == 3) break;
			}

				
		}
		else if ((simpdg == 221) && (particleMap.begin()->second->EndProcess() == "Decay"))
		{	
			return true; 
		}	
		
		return kaonmodedcy;
	}

	return false;
}

bool pid::ProfilePID::prepareEvent(art::Event const & evt)
{
	fRunNum = evt.run(); fEventNum = evt.event();
	fBestView = -1; fTrkIdx = -1;
	fdEdx = -1.0; fRange = -1.0;

	if (!evt.isRealData() && !evselect(evt)) return false;

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
		
			/// choose the index of the interesitng track
			size_t id = 0; int minz = 300;
			for (size_t t = 0; t < hitFromTrk.size(); t++)
			{
				auto const & trk = (*fTrkListHandle)[t];
				if (trk.End().Z() < minz) 
				{
					minz = trk.End().Z(); id = t;
				}
				if (trk.Vertex().Z() < minz)
				{
					minz = trk.Vertex().Z(); id = t;
				}
			}

				double dvtx = 0.0;
				if (hasPfp)
				{
					if (fDirection) // try to get the starting vertex (the way of searching may change...)
					{
						auto const & trk = (*fTrkListHandle)[id];
						auto pfps = pfpFromTrk.at(id);
						if (!pfps.empty())
						{
							auto vtxs = vtxFromPfp.at(pfps.front().key());
							if (!vtxs.empty())
							{
								double vtxpos[3];
								vtxs.front()->XYZ(vtxpos);
								TVector3 vtx3d(vtxpos[0], vtxpos[1], vtxpos[2]);
								dvtx = (vtx3d - trk.Vertex()).Mag();							
							}
						}
					}
					else
					{
						
						auto const & trk = (*fTrkListHandle)[id];
						auto pfps = pfpFromTrk.at(id);
						
						if (!pfps.empty())
						{
							
						}

						if (trk.End().Z() < trk.Vertex().Z()) fFlip = 1;
						
						if ((trk.End().Z() > 1.0) && (trk.End().Z() > 1.0)) return false;
					}
				}

				auto vhit = hitFromTrk.at(id);
				auto vmeta = hitFromTrk.data(id);

				fTrk2VtxDistMap[id] = dvtx;
				
				for (size_t h = 0; h < vhit.size(); ++h)
				{
					size_t view = vhit[h]->WireID().Plane;

					if (view != geo::kZ) continue;					

					size_t idx = vmeta[h]->Index();
					double tdrift = vhit[h]->PeakTime();
					double dx = vmeta[h]->Dx();
					double dqadc = vhit[h]->Integral();
					int wire = vhit[h]->WireID().Wire;

					double dq = fCalorimetryAlg.ElectronsFromADCArea(dqadc, view);
					dq *= fCalorimetryAlg.LifetimeCorrection(tdrift) * (fElectronsToGeV * 1000); // *** note: T0 = 0

					fTrk2InfoMap[id][view].emplace_back(idx, dx, dq, wire);
				}
		}
	}
	if (fTrk2InfoMap.empty()) return false;
	else return true;
}

bool pid::ProfilePID::has(const std::vector<int>& v, int i) const
{
	for (auto c : v) if (c == i) return true;
  return false;
}

void pid::ProfilePID::analyze(art::Event const & evt)
{
	if (!prepareEvent(evt)) return;

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
		}
		else writeTrackInfo(t, dEdx, range); // save training patterns
	}
}

void pid::ProfilePID::make_dEdx(std::vector< double > & dEdx, std::vector< double > & range,
	const std::vector< pid::bHitInfo > & hits, double dvtx, double rmax) const // empty hits is not taken into account.
{
	dEdx.clear(); range.clear();

	int i0, i1, di;
	if (fDirection || fFlip) { i0 = 0; i1 = hits.size(); di = 1; }
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

/*simb::MCParticle const * pid::ProfilePID::getMCParticle(
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
	stopping = true;
	double ek = particle->EndE() - particle->Mass();

	if ((ek < 0.001) || (particle->EndProcess() == "FastScintillation")) stopping = true;
	else stopping = false;

	std::cout << particle->PdgCode() << " ek:" << ek
		<< " process: " << particle->EndProcess()
		<< " clean:" << clean << " stopping:" << stopping
		<< std::endl;

	return particle;
}*/

DEFINE_ART_MODULE(pid::ProfilePID)
