include_directories ( ${PROJECT_SOURCE_DIR} )

set(Calorimtery_HEADERS
     TrackCalorimetryAlg.h
     )

add_library(Calorimetry SHARED
     ${Calorimtery_HEADERS}
     TrackCalorimetryAlg.cxx
)

target_link_libraries(Calorimetry
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
     ${ROOT_BASIC_LIB_LIST}
     )

art_add_module(Calorimetry_module Calorimetry_module.cc)

art_add_module(BezierCalorimetry_module BezierCalorimetry_module.cc)

art_add_module(GeneralCalorimetry_module GeneralCalorimetry_module.cc)

art_add_module(PrintCalorimetry_module PrintCalorimetry_module.cc)

art_add_module(TrackCalorimetry_module TrackCalorimetry_module.cc)

install(TARGETS
     Calorimetry
     Calorimetry_module
     BezierCalorimetry_module
     GeneralCalorimetry_module
     PrintCalorimetry_module
     TrackCalorimetry_module
     EXPORT laranaLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )

install(FILES ${Calorimtery_HEADERS} DESTINATION 
	${CMAKE_INSTALL_INCLUDEDIR}/Calorimetry COMPONENT Development)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)

