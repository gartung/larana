
set(CosmicRemoval_HEADERS
     HitTagAssociatorAlg.h
     BeamFlashTrackMatchTaggerAlg.h
     )


add_library(CosmicRemoval SHARED
     HitTagAssociatorAlg.cxx
     BeamFlashTrackMatchTaggerAlg.cxx
     )

target_link_libraries(CosmicRemoval
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

art_add_module(BeamFlashTrackMatchTagger_module BeamFlashTrackMatchTagger_module.cc)

art_add_module(CRHitRemoval_module CRHitRemoval_module.cc)

art_add_module(CosmicClusterTagger_module CosmicClusterTagger_module.cc)

art_add_module(CosmicPFParticleTagger_module CosmicPFParticleTagger_module.cc)

art_add_module(CosmicRemovalAna_module CosmicRemovalAna_module.cc)

art_add_module(CosmicTrackTagger_module CosmicTrackTagger_module.cc)


install(TARGETS
     CosmicRemoval
     BeamFlashTrackMatchTagger_module
     CRHitRemoval_module
     CosmicClusterTagger_module
     CosmicPFParticleTagger_module
     CosmicRemovalAna_module
     CosmicTrackTagger_module
     EXPORT laranaLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )

install(FILES ${CosmicRemoval_HEADERS} DESTINATION
	${CMAKE_INSTALL_INCLUDEDIR}/CosmicRemoval COMPONENT Development)	

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)
