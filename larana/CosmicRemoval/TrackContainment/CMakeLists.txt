# basic source code CMakeLists.txt

art_make( 
  LIB_LIBRARIES
  larcorealg_Geometry
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  ${FHICLCPP}
  ${ROOT_LIBRARIES}
  cetlib cetlib_except
  ${BOOST_SYSTEM_LIBRARY} #why the hell do I need this?
  )



install_headers()
install_source()
