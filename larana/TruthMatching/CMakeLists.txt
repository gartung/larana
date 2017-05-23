include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( 
	  
	  LIB_LIBRARIES
            larsim_MCCheater_BackTracker_service
			nusimdata_SimulationBase
                        lardataobj_AnalysisBase
			larreco_Calorimetry
                        lardata_Utilities
                        larcore_Geometry
            lardataobj_RecoBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        cetlib cetlib_except
		        ${CLHEP}
		        ${ROOT_BASIC_LIB_LIST}
                        ${Boost_SYSTEM_LIBRARY}
         )

install_headers()
install_fhicl()
install_source()
