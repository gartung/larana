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
 * track direction at that point less than maxAngleRad.
 *
 * distance is set to the distance between the track endpoint
 * and the mcParticle start point.
 */
const art::Ptr<recob::Track> 
mctrue::StartPosDirMatcherAlg::getBestMatch(
            simb::MCParticle const& mcParticle, 
            std::vector<art::Ptr<recob::Track>> const& recoTracks,
            double const& maxAngleRad, double& distance)
{
  const TVector3 mcpStartPos = mcParticle.Position().Vect();
  const TVector3 mcpStartMom = mcParticle.Momentum().Vect();

  auto bestMatch = recoTracks.end();
  double bestDistance = 1e9; 
  for(auto track = recoTracks.begin(); track != recoTracks.end(); track++)
  {
    // First from the front of the track
    const TVector3 vertex = (*track)->Vertex();
    double vertexAngle = mcpStartMom.Angle(vertex);
    if (vertexAngle > TMath.PiOver2())
    {
      vertexAngle -= TMath.PiOver2();
    }
    if (vertexAngle <= maxAngleRad)
    {
      const double vertexDistance = (mcpStartPos-vertex).Mag();
      if (vertexDistance < bestDistance)
      {
        bestDistance = vertexDistance;
        bestMatch = track;
      }
    }
    // Then end of track
    const TVector3 trkEnd = (*track)->End();
    const double endAngle = mcpStartMom.Angle(trkEnd);
    if (endAngle > TMath.PiOver2())
    {
      endAngle -= TMath.PiOver2();
    }
    if (endAngle <= maxAngleRad)
    {
      const double endDistance = (mcpStartPos-trkEnd).Mag();
      if (endDistance < bestDistance)
      {
        bestDistance = endDistance;
        bestMatch = track;
      }
    }
  } //loop over recob::Tracks

  if(bestMatch != recoVect.end())
  {
    return *bestMatch; 
  }
  return art::Ptr<recob::Track>(); //return an invalid art::Ptr
}
