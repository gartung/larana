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
            bool extrapolate)
{
  const size_t ntrajectory = trajectory.size();
  const int iLastTrajectory = ntrajectory-1;
  if (ntrajectory == 0)
  {
    distance = -1.;
    return TVector3();
  }
  else if (ntrajectory == 1)
  {
    distance = (point-trajectory.Position(0).Vect()).Mag();
    interpolatedMomentum = trajectory.Momentum(0);
    return trajectory.Position(0).Vect();
  }
  const int iClosest = findClosestTrajPoint(trajectory,point);
  if (iClosest < 0)
  {
    distance = -1.;
    return TVector3();
  }
  else if (iClosest == 0)
  {
    return interpExtrapAtBeginning(trajectory,point,distance,interpolatedMomentum,extrapolate);
  }
  else if (iClosest == iLastTrajectory)
  {
    return interpExtrapAtEnd(trajectory,point,distance,interpolatedMomentum,extrapolate);
  }
  else
  {
    return interplateInMiddle(trajectory,point,distance,interpolatedMomentum,iClosest);
  }
}

/// Find point of closest approach when closest traj point in middle of traj points
/**
 * requires that the closest trajectory point has already been found, iClosest
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::interplateInMiddle(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            const int iClosest)
{
  double locationMeasureBefore = 0.;
  double locationMeasureAfter = 0.;
  TVector3 segmentBefore = trajectory.Position(iClosest-1).Vect();
  TVector3 segmentCenter = trajectory.Position(iClosest).Vect();
  TVector3 segmentAfter = trajectory.Position(iClosest+1).Vect();
  TVector3 closestBefore = segmentPointOfClosestApproach(segmentBefore,segmentCenter,point,locationMeasureBefore);
  TVector3 closestAfter = segmentPointOfClosestApproach(segmentAfter,segmentCenter,point,locationMeasureAfter);
  if (locationMeasureBefore < 0)
  {
    throw cet::exception("IvalidValue","locationMeasureBefore < 0 when segmentCenter should be closer");
  }
  if (locationMeasureBefore > 1)
  {
      closestBefore = segmentCenter;
  }
  if (locationMeasureAfter < 0)
  {
      closestAfter = segmentCenter;
  }
  if (locationMeasureAfter > 1)
  {
    throw cet::exception("IvalidValue","locationMeasureAfter > 1 when segmentCenter should be closer");
  }
  double distanceBefore = (point-closestBefore).Mag();
  double distanceAfter = (point-closestAfter).Mag();
  if (distanceBefore < distanceAfter)
  {
    if(locationMeasureBefore > 1)
    {
      interpolatedMomentum = trajectory.Momentum(iClosest);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                              closestBefore,
                              trajectory.Position(iClosest-1).Vect(),
                              trajectory.Position(iClosest).Vect(),
                              trajectory.Momentum(iClosest-1),
                              trajectory.Momentum(iClosest));
    }
    distance = distanceBefore;
    return closestBefore;
  }
  else
  {
    if(locationMeasureAfter < 0)
    {
      interpolatedMomentum = trajectory.Momentum(iClosest);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                              closestAfter,
                              trajectory.Position(iClosest).Vect(),
                              trajectory.Position(iClosest+1).Vect(),
                              trajectory.Momentum(iClosest),
                              trajectory.Momentum(iClosest+1));
    }
    distance = distanceAfter;
    return closestAfter;
  }
}

/// Find point of closest approach when closest traj point is the first
/**
 * Assumes the closest trajectory point is the first
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::interpExtrapAtBeginning(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            bool extrapolate)
{
  double locationMeasure = 0.;
  TVector3 segment1 = trajectory.Position(0).Vect();
  TVector3 segment2 = trajectory.Position(1).Vect();
  TVector3 thisPoint = segmentPointOfClosestApproach(segment1,segment2,point,locationMeasure);
  if ((!extrapolate) && locationMeasure < 0.)
  {
    distance = (point-trajectory.Position(0).Vect()).Mag();
    interpolatedMomentum = trajectory.Momentum(0);
    return trajectory.Position(0).Vect();
  }
  else if (locationMeasure > 1.)
  {
    throw cet::exception("IvalidValue","locationMeasure > 1 when segment1 should be closer");
  }
  else
  {
    distance = (point-thisPoint).Mag();
    if (extrapolate && locationMeasure < 0.)
    {
      interpolatedMomentum = trajectory.Momentum(0);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                                thisPoint,
                                trajectory.Position(0).Vect(),
                                trajectory.Position(1).Vect(),
                                trajectory.Momentum(0),
                                trajectory.Momentum(1));
    }
    return thisPoint;
  }
}

/// Find point of closest approach when closest traj point is the last
/**
 * Assumes the closest trajectory point is the last
 */
const TVector3 
mctrue::TrajectoryInterpExtrapAlg::interpExtrapAtEnd(
            const simb::MCTrajectory& trajectory,
            const TVector3& point,
            double& distance,
            TLorentzVector& interpolatedMomentum,
            bool extrapolate)
{
  double locationMeasure = 0.;
  const size_t ntrajectory = trajectory.size();
  const int iLastTrajectory = ntrajectory-1;
  TVector3 segment1 = trajectory.Position(iLastTrajectory).Vect();
  TVector3 segment2 = trajectory.Position(iLastTrajectory-1).Vect();
  TVector3 thisPoint = segmentPointOfClosestApproach(segment1,segment2,point,locationMeasure);
  if ((!extrapolate) && locationMeasure < 0.)
  {
    distance = (point-trajectory.Position(iLastTrajectory).Vect()).Mag();
    interpolatedMomentum = trajectory.Momentum(iLastTrajectory);
    return trajectory.Position(iLastTrajectory).Vect();
  }
  else if (locationMeasure > 1.)
  {
    throw cet::exception("IvalidValue","locationMeasure > 1 when segment1 should be closer");
  }
  else
  {
    distance = (point-thisPoint).Mag();
    if (extrapolate && locationMeasure < 0.)
    {
      interpolatedMomentum = trajectory.Momentum(iLastTrajectory);
    }
    else
    {
      interpolatedMomentum = interpolateMomentum(
                              thisPoint,
                              trajectory.Position(iLastTrajectory).Vect(),
                              trajectory.Position(iLastTrajectory-1).Vect(),
                              trajectory.Momentum(iLastTrajectory),
                              trajectory.Momentum(iLastTrajectory-1));
    }
    return thisPoint;
  }
}

/// Finds the trajectory point that is closest
/**
 * Returns -1 if none found, and the index of the closest point otherwise
 */
int
mctrue::TrajectoryInterpExtrapAlg::findClosestTrajPoint(
            const simb::MCTrajectory& trajectory,
            const TVector3& point)
{
  const size_t ntrajectory = trajectory.size();
  int result = -1;
  double closestDistance = 1e15;
  for(size_t iPoint=0; iPoint<ntrajectory; iPoint++)
  {
    double thisDistance = (point-trajectory.Position(iPoint).Vect()).Mag();
    if (thisDistance < closestDistance)
    {
      result = iPoint;
      closestDistance = thisDistance;
    }
  }
  return result;
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
            double& locationMeasure)
{
   const TVector3 u = point-segment1;
   const TVector3 v = segment2-segment1;
   const double uMag = u.Mag();
   locationMeasure = u.Dot(v)/uMag;
   const double vMagInv = 1./v.Mag();
   const TVector3 vhat = v * vMagInv;
   const TVector3 result = segment1 + u.Dot(vhat)*vhat;
   return result;
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
    const double pointDistance = (point2-point1).Mag();
    const double approachDistance = (pointOfClosestApproach-point1).Mag();
    const double distanceFraction = pointDistance/approachDistance;
    const TLorentzVector momDiff = momentum2-momentum1;
    return distanceFraction*momDiff + momentum1;
}

