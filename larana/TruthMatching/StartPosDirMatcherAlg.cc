/*!
 * Title:   Start Position & Direction Matcher Algorithim Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for matching MCParticles to reconstructed tracks 
 *              based on the start position and direction of the track and 
 *              MCParticle.
 * Input:       simb::MCParticle, vector<art::Ptr<recob::Track>>, and 
 *              maximum angle for match (double)
 * Output:      Matched track (art::Ptr<recob::Track>) and angle in 
 *              radians between start direction of track and MCParticle 
 *              (double)
*/

#include "StartPosDirMatcherAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TMath.h"
#include <algorithm>

/// Find best match track given MCParticle and set of tracks
/**
 * Find the best matching track to the given mcParticle.
 * The best match is found by finding the recoTrack with
 * one end closest to the start of the mcParticle trajectory
 * and with angle between the trajectory direction and
 * track direction at that point less than maxAngleDeg.
 * maxAngleDeg is the minimum angle between the two lines
 * i.e. always <= 90 deg.
 *
 * distance is set to the distance between the track endpoint
 * and the mcParticle start point.
 *
 * angleDeg is set to the angle between the track direction at
 * the matched endpoint and the mcParticle start direction.
 * It is always <= 90 deg.
 */
const art::Ptr<recob::Track> 
mctrue::StartPosDirMatcherAlg::getBestMatch(
            simb::MCParticle const& mcParticle, 
            std::vector<art::Ptr<recob::Track>> const& recoTracks,
            double const& maxAngleDeg, double& distance, double& angleDeg)
{
  distance = 1e9;
  if(mcParticle.NumberTrajectoryPoints() == 0)
  {
    return art::Ptr<recob::Track>(); //return an invalid art::Ptr
  }
  const TVector3 mcpStartPos = mcParticle.Position().Vect();
  const TVector3 mcpStartMom = mcParticle.Momentum().Vect();
  const double maxAngleRad = maxAngleDeg*TMath::Pi()/180.;

  auto bestMatch = recoTracks.end();
  for(auto track = recoTracks.begin(); track != recoTracks.end(); track++)
  {
    // First from the front of the track
    const TVector3 vertex = (*track)->Vertex();
    double vertexAngle = mcpStartMom.Angle(vertex);
    if (vertexAngle > TMath::PiOver2())
    {
      vertexAngle -= TMath::PiOver2();
    }
    if (vertexAngle <= maxAngleRad)
    {
      const double vertexDistance = (mcpStartPos-vertex).Mag();
      if (vertexDistance < distance)
      {
        distance = vertexDistance;
        angleDeg = vertexAngle * 180 / TMath::Pi();
        bestMatch = track;
      }
    }
    // Then end of track
    const TVector3 trkEnd = (*track)->End();
    double endAngle = mcpStartMom.Angle(trkEnd);
    if (endAngle > TMath::PiOver2())
    {
      endAngle -= TMath::PiOver2();
    }
    if (endAngle <= maxAngleRad)
    {
      const double endDistance = (mcpStartPos-trkEnd).Mag();
      if (endDistance < distance)
      {
        distance = endDistance;
        angleDeg = vertexAngle * 180 / TMath::Pi();
        bestMatch = track;
      }
    }
  } //loop over recob::Tracks

  if(bestMatch != recoTracks.end())
  {
    return *bestMatch; 
  }
  return art::Ptr<recob::Track>(); //return an invalid art::Ptr
}
