#define BOOST_TEST_MODULE ( TrajectoryInterpExtrapAlg_test )
#include "cetlib/quiet_unit_test.hpp"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "larana/TruthMatching/TrajectoryInterpExtrapAlg.h"


#define checkTVec3Close(a, b, tolerance) { \
  BOOST_REQUIRE_CLOSE(a.X(),b.X(),tolerance); \
  BOOST_REQUIRE_CLOSE(a.Y(),b.Y(),tolerance); \
  BOOST_REQUIRE_CLOSE(a.Z(),b.Z(),tolerance); \
}
//void checkTVec3Close(const TVector3 & a, const TVector3 & b, double tolerance)
//{
//  BOOST_REQUIRE_CLOSE(a.X(),b.X(),tolerance);
//  BOOST_REQUIRE_CLOSE(a.Y(),b.Y(),tolerance);
//  BOOST_REQUIRE_CLOSE(a.Z(),b.Z(),tolerance);
//}

#define checkTLVecClose(a,b,tolerance) { \
  BOOST_REQUIRE_CLOSE(a.X(),b.X(),tolerance); \
  BOOST_REQUIRE_CLOSE(a.Y(),b.Y(),tolerance); \
  BOOST_REQUIRE_CLOSE(a.Z(),b.Z(),tolerance); \
  BOOST_REQUIRE_CLOSE(a.T(),b.T(),tolerance); \
}

//void checkTLVecClose(const TLorentzVector & a, const TLorentzVector & b, double tolerance)
//{
//  BOOST_REQUIRE_CLOSE(a.X(),b.X(),tolerance);
//  BOOST_REQUIRE_CLOSE(a.Y(),b.Y(),tolerance);
//  BOOST_REQUIRE_CLOSE(a.Z(),b.Z(),tolerance);
//  BOOST_REQUIRE_CLOSE(a.T(),b.T(),tolerance);
//}

struct TrajectoryInterpExtrapAlgFixture{

  TrajectoryInterpExtrapAlgFixture() : myTrajAlg() 
  {

    //mf::StartMessageFacility( mf::MessageFacilityService::SingleThread,
    //                             mf::MessageFacilityService::logConsole());
  };
  mctrue::TrajectoryInterpExtrapAlg myTrajAlg;

};

double const tolerance = 1e-6;

BOOST_FIXTURE_TEST_SUITE(TrajectoryInterpExtrapAlg_test, TrajectoryInterpExtrapAlgFixture)

BOOST_AUTO_TEST_CASE(checkNoTrajPoints)
{
  simb::MCTrajectory mcTraj;
  TLorentzVector interpMom;
  TVector3 point;
  double distance;
  TVector3 result;

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance < 0);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance < 0);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance < 0);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance < 0);
  
}

BOOST_AUTO_TEST_CASE(checkOneTrajPoints)
{
  simb::MCTrajectory mcTraj;
  TLorentzVector firstTrajPoint(0,0,0,0);
  TLorentzVector firstTrajPointMom(1,0,0,2);
  mcTraj.push_back(firstTrajPoint,firstTrajPointMom);
  TLorentzVector interpMom;
  TVector3 point;
  double distance;
  TVector3 result;

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,firstTrajPointMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,firstTrajPointMom,1e-3);

  point = TVector3(1.,0.,0.);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,1.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,1.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,firstTrajPointMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,1.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,1.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,firstTrajPointMom,1e-3);
  
}

BOOST_AUTO_TEST_CASE(checkTwoTrajPoints)
{
  simb::MCTrajectory mcTraj;
  TLorentzVector firstTrajPoint(0,0,0,0);
  TLorentzVector firstTrajPointMom(1,0,0,5);
  TLorentzVector secondTrajPoint(1,0,0,1);
  TLorentzVector secondTrajPointMom(0,0,0,1);
  mcTraj.push_back(firstTrajPoint,firstTrajPointMom);
  mcTraj.push_back(secondTrajPoint,secondTrajPointMom);
  TLorentzVector interpMom;
  TVector3 point;
  double distance;
  TVector3 result;

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,firstTrajPointMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,firstTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,firstTrajPointMom,1e-3);

  point = TVector3(1.,0.,0.);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,secondTrajPoint.Vect(),1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,secondTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,secondTrajPointMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,secondTrajPoint.Vect(),1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,0.,1e-3);
  checkTVec3Close(result,secondTrajPoint.Vect(),1e-3);
  checkTLVecClose(interpMom,secondTrajPointMom,1e-3);

  point = TVector3(0.75,5.,0.);
  TVector3 correctPos(0.75,0,0);
  TLorentzVector correctMom(0.25,0,0,2);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  //mf::LogDebug ldb("interpolateMomentum args: ");
  // ldb << "result: " << result.X()
  //     << ", " << result.Y() << ", " << result.Z() << " ";
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,5.,1e-3);
  checkTVec3Close(result,correctPos,1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,5.,1e-3);
  checkTVec3Close(result,correctPos,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,5.,1e-3);
  checkTVec3Close(result,correctPos,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,5.,1e-3);
  checkTVec3Close(result,correctPos,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);
  
  point = TVector3(-1.,5.,0.);
  TVector3 correctPosInterp = firstTrajPoint.Vect();
  TVector3 correctPosExtrap(-1,0.,0);
  correctMom = firstTrajPointMom;
  double correctDistanceInterp = (point-correctPosInterp).Mag();
  double correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  //mf::LogDebug ldb("interpolateMomentum args: ");
  // ldb << "result: " << result.X()
  //     << ", " << result.Y() << ", " << result.Z() << " ";
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  point = TVector3(6.,5.,0.125);
  correctPosInterp = secondTrajPoint.Vect();
  correctPosExtrap = TVector3(6.,0.,0.);
  correctMom = secondTrajPointMom;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  
  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,true);
  BOOST_CHECK(distance >= 0.);
  //mf::LogDebug ldb("interpolateMomentum args: ");
  // ldb << "result: " << result.X()
  //     << ", " << result.Y() << ", " << result.Z() << " ";
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);
  
}

BOOST_AUTO_TEST_CASE(checkLongLineTraj)
{
  simb::MCTrajectory mcTraj;
  size_t nPoints = 1000;
  for (size_t i=0;i<nPoints;i++)
  {
    mcTraj.push_back(TLorentzVector(i,0,0,0),TLorentzVector(i,0,0,i));
  }
  TVector3 point;
  TLorentzVector interpMom;
  double distance;
  TVector3 result;

  /////////

  TVector3 correctPosInterp = mcTraj.Position(0).Vect();
  TVector3 correctPosExtrap = correctPosInterp;
  TLorentzVector correctMom = mcTraj.Momentum(0);
  double correctDistanceInterp = (point-correctPosInterp).Mag();
  double correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  /////////

  point = mcTraj.Position(nPoints-1).Vect()+TVector3(0,5,10);
  correctPosInterp = mcTraj.Position(nPoints-1).Vect();
  correctPosExtrap = correctPosInterp;
  correctMom = mcTraj.Momentum(nPoints-1);
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);
  
  /////////

  point = mcTraj.Position((nPoints-1)/2).Vect()+TVector3(0,5,10);
  correctPosInterp = mcTraj.Position((nPoints-1)/2).Vect();
  correctPosExtrap = correctPosInterp;
  correctMom = mcTraj.Momentum((nPoints-1)/2);
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  /////////

  point = mcTraj.Position(0).Vect()+TVector3(-10,5,10);
  correctPosInterp = mcTraj.Position(0).Vect();
  correctPosExtrap =  mcTraj.Position(0).Vect()+TVector3(-10,0,0);
  correctMom = mcTraj.Momentum(0);
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);
  
  /////////

  point = mcTraj.Position(nPoints-1).Vect()+TVector3(100,-500,10);
  correctPosInterp = mcTraj.Position(nPoints-1).Vect();
  correctPosExtrap =  mcTraj.Position(nPoints-1).Vect()+TVector3(100,0,0);
  correctMom = mcTraj.Momentum(nPoints-1);
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);
  
  /////////

  point = mcTraj.Position((nPoints-1)/2).Vect()+TVector3(0.2,5,10);
  correctPosInterp = mcTraj.Position((nPoints-1)/2).Vect()+TVector3(0.2,0,0);
  correctPosExtrap = correctPosInterp;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);

  /////////

  point = mcTraj.Position((nPoints-1)/2).Vect()+TVector3(-0.2,5,10);
  correctPosInterp = mcTraj.Position((nPoints-1)/2).Vect()+TVector3(-0.2,0,0);
  correctPosExtrap = correctPosInterp;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);

  /////////

  point = mcTraj.Position(0).Vect()+TVector3(0.45,5,10);
  correctPosInterp = mcTraj.Position(0).Vect()+TVector3(0.45,0,0);
  correctPosExtrap = correctPosInterp;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  
  /////////

  point = mcTraj.Position(nPoints-1).Vect()+TVector3(-0.49,-500,10);
  correctPosInterp = mcTraj.Position(nPoints-1).Vect()+TVector3(-0.49,0,0);
  correctPosExtrap = correctPosInterp;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  
}

BOOST_AUTO_TEST_CASE(checkSharpKinkTraj)
{
  simb::MCTrajectory mcTraj;
  TLorentzVector firstTrajPoint(0,0,0,0);
  TLorentzVector firstTrajPointMom(1,0,0,1);
  TLorentzVector secondTrajPoint(1,1,0,1);
  TLorentzVector secondTrajPointMom(2,0,0,2);
  TLorentzVector thirdTrajPoint(2,0,0,2);
  TLorentzVector thirdTrajPointMom(3,0,0,3);
  mcTraj.push_back(firstTrajPoint,firstTrajPointMom);
  mcTraj.push_back(secondTrajPoint,secondTrajPointMom);
  mcTraj.push_back(thirdTrajPoint,thirdTrajPointMom);
  TVector3 point;
  TLorentzVector interpMom;
  double distance;
  TVector3 result;

  /////////

  point = TVector3(1,2,0);
  TVector3 correctPosInterp = secondTrajPoint.Vect();
  TVector3 correctPosExtrap = correctPosInterp;
  TLorentzVector correctMom = secondTrajPointMom;
  double correctDistanceInterp = (point-correctPosInterp).Mag();
  double correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  /////////

  point = TVector3(1.1,2,0);
  correctPosInterp = secondTrajPoint.Vect();
  correctPosExtrap = correctPosInterp;
  correctMom = secondTrajPointMom;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  /////////

  point = TVector3(0.9,2,0);
  correctPosInterp = secondTrajPoint.Vect();
  correctPosExtrap = correctPosInterp;
  correctMom = secondTrajPointMom;
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  /////////

  point = TVector3(0,1.,0);
  correctPosInterp = TVector3(0.5,0.5,0.);
  correctPosExtrap = correctPosInterp;
  correctMom = TLorentzVector(1.5,0,0,1.5);
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  /////////

  point = TVector3(2.,1.,0);
  correctPosInterp = TVector3(1.5,0.5,0.);
  correctPosExtrap = correctPosInterp;
  correctMom = TLorentzVector(2.5,0,0,2.5);
  correctDistanceInterp = (point-correctPosInterp).Mag();
  correctDistanceExtrap = (point-correctPosExtrap).Mag();

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceInterp,1e-3);
  checkTVec3Close(result,correctPosInterp,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

  result = myTrajAlg.pointOfClosestApproach(mcTraj,point,distance,interpMom,true);
  BOOST_CHECK(distance >= 0.);
  BOOST_REQUIRE_CLOSE(distance,correctDistanceExtrap,1e-3);
  checkTVec3Close(result,correctPosExtrap,1e-3);
  checkTLVecClose(interpMom,correctMom,1e-3);

}


BOOST_AUTO_TEST_SUITE_END()
