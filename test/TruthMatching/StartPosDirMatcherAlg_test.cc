#define BOOST_TEST_MODULE ( StartPosDirMatcherAlg_test )
#include "cetlib/quiet_unit_test.hpp"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "larana/TruthMatching/StartPosDirMatcherAlg.h"

struct StartPosDirMatcherAlgFixture{

  StartPosDirMatcherAlgFixture() : matcher() 
  {

    //mf::StartMessageFacility( mf::MessageFacilityService::SingleThread,
    //                             mf::MessageFacilityService::logConsole());
  };
  mctrue::StartPosDirMatcherAlg matcher;

};

double const tolerance = 1e-6;

BOOST_FIXTURE_TEST_SUITE(StartPosDirMatcherAlg_test, StartPosDirMatcherAlgFixture)

BOOST_AUTO_TEST_CASE(checkEmpty)
{
  simb::MCParticle mcPart;
  std::vector<art::Ptr<recob::Track>> tracks;
  double maxAngleRad=.1;
  double distance;
  double angle;
  art::Ptr<recob::Track> result;

  result = matcher.getBestMatch(mcPart,tracks,maxAngleRad);
  BOOST_CHECK(distance >= 1e9);
  BOOST_CHECK_EQUAL(result,art::Ptr<recob::Track>());
  result = matcher.getBestMatch(mcPart,tracks,maxAngleRad,distance,angle);
  BOOST_CHECK(distance >= 1e9);
  BOOST_CHECK_EQUAL(result,art::Ptr<recob::Track>());
}

BOOST_AUTO_TEST_SUITE_END()
