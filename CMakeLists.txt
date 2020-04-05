# ======================================================================
#  larana main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larana/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 3.11 or later
cmake_minimum_required (VERSION 3.11)

project(larana VERSION 08.10.03 LANGUAGES CXX)

set(fcl_dir job)
set(fw_dir config_data)
set(gdml_dir gdml)

find_package(cetmodules REQUIRED)
include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package( larreco REQUIRED )


# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(larana)

# tests
add_subdirectory(test)
add_subdirectory(cmakeConfig)
# ups - table and config files
#add_subdirectory(ups)

# packaging utility
#include(UseCPack)

