/*!
 * Title:   Trjectory Interpolation and Extrapolation Algorithms Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for interpolating and extrapolating trajectories
 * Input:       simb::MCTrajectory trajectory, TVector3 point, bool extraploate or not
 * Output:      TVector3 point of closest approach, distance of closest approach, and interplated momentum
*/

#include "TrajectoryInterpExtrapAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

mctrue::TrajectoryInterpExtrapAlg::BackTrackMatcherAlg(fhicl::ParameterSet const& pset)
{
  reconfigure(pset);
  return;
} //end constructor

void 
mctrue::TrajectoryInterpExtrapAlg::reconfigure(fhicl::ParameterSet const& pset)
{
  return;
}

/// Finds the point of closest approach of the trajectory to the point
/**
 * distance is set to the distance between the point of closest
 * approach and the input point. If extrapolate is false and
 * distance is large, this may mean that the point is nowhere near
 * the trajectory.
 *
 * If extrapolate is true, will try projecting the trajectory
 * start point backwards and end point forwards in search of
 * a point of closer approach.
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::pointOfClosestApproach(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            bool extrapolate=false);
{
  return TVector3();
}

/// Finds the point of closest approach of the trajectory to the point
/**
 * distance is set to the distance between the point of closest
 * approach and the input point. If extrapolate is false and
 * distance is large, this may mean that the point is nowhere near
 * the trajectory.
 *
 * interpolatedMomentum will be set to the approximate four-momentum 
 * at the point of closest approach. If the point of closest approach
 * is between trajectory points, the momentum is computed by linearly
 * interpolating the momentum between them. Otherwise, it is set to
 * the momentum of the nearest trajectory point.
 *
 * If extrapolate is true, will try projecting the trajectory
 * start point backwards and end point forwards in search of
 * a point of closer approach.
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::pointOfClosestApproach(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            bool extrapolate=false);
{
  return TVector3();
}

/// Finds the point of closest approach to the line going through the segment
/**
 * The resulting vector is extrapolated out along the line
 * going through the points
 *
 * The values of locationMeasure are set such that:
 *    location Measure < 0: segment1 is the closest point on the segment
 *    location Measure > 1: segment2 is the closest point on the segment
 *    otherwise: the returned point is the closest point on the segment
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::segmentPointOfClosestApproach(
            const TVector3& segment1,
            const TVector3& segment2,
            const TVector3& point,
            double& locationMeasure);
{
  return TVector3();
}

/// Finds the point of closest approach, no extrapolation
/**
 * distance is set to the distance between the point of closest
 * approach and the input point. If extrapolate is false and
 * distance is large, this may mean that the point is nowhere near
 * the trajectory.
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::pointOfClosestApproachNoExtrap(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            TLorentzVector& interpolatedMomentum,
            double& distance);
{
  return TVector3();
}

/// Finds the point of closest approach, with extrapolation
/**
 * distance is set to the distance between the point of closest
 * approach and the input point. If extrapolate is false and
 * distance is large, this may mean that the point is nowhere near
 * the trajectory.
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::pointOfClosestApproachDoExtrap(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            TLorentzVector& interpolatedMomentum,
            double& distance);
{
  return TVector3();
}

/// Find the momentum at the point of closest approach
/**
 * distance is set to the distance between the point of closest
 * approach and the input point. If extrapolate is false and
 * distance is large, this may mean that the point is nowhere near
 * the trajectory.
 */
const TLorentzVector 
mctrue::TrajectoryInterpExtrapAlg::interpolateMomentum(
            const TVector3& pointOfClosestApproach,
            const TVector3& point1,
            const TVector3& point2,
            const TLorentzVector& momentum1,
            const TLorentzVector& momentum2)
{
  return TLorentzVector();
}

