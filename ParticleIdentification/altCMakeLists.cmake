include_directories ( ${PROJECT_SOURCE_DIR} )

set(ParticleIdentification_HEADERS
     PIDAAlg.h
     )

add_library(ParticleIdentification SHARED
     ${ParticleIdentification_HEADERS}
     PIDAAlg.cxx
     )

target_link_libraries(ParticleIdentification
     larsoft::AnalysisBase
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

art_add_module( Chi2ParticleID_module Chi2ParticleID_module.cc )

art_add_module( PIDAAnalyzer_module PIDAAnalyzer_module.cc )

install(TARGETS
     ParticleIdentification
     Chi2ParticleID_module
     PIDAAnalyzer_module
     EXPORT laranaLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )

install(FILES ${ParticleIdentification_HEADERS} DESTINATION 
	${CMAKE_INSTALL_INCLUDEDIR}/ParticleIdentification COMPONENT Development)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)

