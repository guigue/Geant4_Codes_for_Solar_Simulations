# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(solar_test)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(solar solar.cc ${sources} ${headers})
target_link_libraries(solar ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build solar_test. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SOLAR_SCRIPTS
  solar.in
  solar.out
  init_vis.mac
  run0.mac
  runTest.mac
  runTestEnergy.mac
  vis.mac
  test_solar1.mac
  Solar_field.in
  Solar_field.out
  run25015.mac
  run25035.mac
  run25080.mac
  run115.mac
  run135.mac
  run180.mac
  erun10015.mac
  erun10035.mac
  erun10080.mac
  run515.mac
  run535.mac
  run580.mac
  run50015.mac
  run50035.mac
  run50080.mac
  erun115.mac
  erun135.mac
  erun180.mac
  erun1015.mac
  erun1035.mac
  erun1080.mac
  erun50015.mac
  erun50035.mac
  erun50080.mac
  )

foreach(_script ${SOLAR_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
#
add_custom_target(solar_test DEPENDS solar)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS solar DESTINATION bin)
