set(OpticalDetector_HEADERS
     AlgoFixedWindow.h
     AlgoPedestal.h
     AlgoThreshold.h
     DefaultOpDetResponse.h
     FlashHypothesis.h
     FlashHypothesisAnaAlg.h
     FlashHypothesisCalculator.h
     FlashHypothesisComparison.h
     FlashHypothesisCreator.h
     FlashUtilities.h
     MicrobooneOpDetResponse.h
     OpDetResponseInterface.h
     OpDigiProperties.h
     OpFlashAlg.h
     OpticalRecoAna.h
     PMTPulseRecoBase.h
     PulseRecoManager.h
     SimPhotonCounter.h
     SimPhotonCounterAlg.h
     )

add_library(OpticalDetector SHARED
     ${OpticalDetector_HEADERS}
     AlgoFixedWindow.cxx
     AlgoPedestal.cxx
     AlgoThreshold.cxx
     FlashHypothesis.cxx
     FlashHypothesisAnaAlg.cxx
     FlashHypothesisCalculator.cxx
     FlashHypothesisComparison.cxx
     FlashHypothesisCreator.cxx
     FlashUtilities.cxx
     OpFlashAlg.cxx
     PMTPulseRecoBase.cxx
     PulseRecoManager.cxx
     SimPhotonCounter.cxx
     SimPhotonCounterAlg.cxx
     )

target_link_libraries(OpticalDetector
     larsoft::RecoBase
     larsoft::Geometry
     larsoft::RecoAlg
     larsoft::Utilities
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
     )

art_add_service(OpDigiProperties_service OpDigiProperties_service.cc)

art_add_service(DefaultOpDetResponse_service DefaultOpDetResponse_service.cc)

art_add_service(MicrobooneOpDetResponse_service MicrobooneOpDetResponse_service.cc)

art_add_module( BeamFlashCompatibilityCheck_module BeamFlashCompatibilityCheck_module.cc )

art_add_module( BoDataFrameInput_module BoDataFrameInput_module.cc )

art_add_module( FIFOHistogramAna_module FIFOHistogramAna_module.cc )

art_add_module( FlashClusterMatch_module FlashClusterMatch_module.cc )

art_add_module( FlashHypothesisAna_module FlashHypothesisAna_module.cc )

art_add_module( FlashPurityCheckAna_module FlashPurityCheckAna_module.cc)

art_add_module( LEDCalibrationAna_module LEDCalibrationAna_module.cc )

art_add_module( OpDigiAna_module OpDigiAna_module.cc )

art_add_module( OpFlashAna_module OpFlashAna_module.cc )

art_add_module( OpFlashFinder_module OpFlashFinder_module.cc )

art_add_module( OpFlashMCTruthAna_module OpFlashMCTruthAna_module.cc )

art_add_module( OpHitAna_module OpHitAna_module.cc )

art_add_module( OpMCDigi_module OpMCDigi_module.cc )

art_add_module( OptDetDigitizer_module OptDetDigitizer_module.cc )

art_add_module( OpticalRecoAna_module OpticalRecoAna_module.cc )

art_add_module( PMTAna_module PMTAna_module.cc )

art_add_module( SimPhotonCounter_module SimPhotonCounter_module.cc )

art_add_module( TrackTimeAssocAna_module TrackTimeAssocAna_module.cc )

art_add_module( TrackTimeAssoc_module TrackTimeAssoc_module.cc )


install(TARGETS
     OpticalDetector
     OpDigiProperties_service
     DefaultOpDetResponse_service
     MicrobooneOpDetResponse_service
     BeamFlashCompatibilityCheck_module
     BoDataFrameInput_module
     FIFOHistogramAna_module
     FlashClusterMatch_module
     FlashHypothesisAna_module
     FlashPurityCheckAna_module
     LEDCalibrationAna_module
     OpDigiAna_module
     OpFlashAna_module
     OpFlashFinder_module
     OpFlashMCTruthAna_module
     OpHitAna_module
     OpMCDigi_module
     OptDetDigitizer_module
     OpticalRecoAna_module
     PMTAna_module
     SimPhotonCounter_module
     TrackTimeAssocAna_module
     TrackTimeAssoc_module
     EXPORT laranaLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )

install(FILES ${OpticalDetector_HEADERS} DESTINATION 
	${CMAKE_INSTALL_INCLUDEDIR}/OpticalDetector COMPONENT Development)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)

install(FILES 
     toyWaveform.txt
     DESTINATION config_data/Optical_Detector 
     COMPONENT Runtime
     )
