# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_test(TrajectoryInterpExtrapAlg_test USE_BOOST_UNIT
			 LIBRARIES larana_TruthMatching
)

