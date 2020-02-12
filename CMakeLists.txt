cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

# If ccache is installed use it
find_program(CCACHE_EXECUTABLE ccache)
if(CCACHE_EXECUTABLE)
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_EXECUTABLE})
	set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_EXECUTABLE})
endif()

project(tpc-rs LANGUAGES CXX)

# By default build shared libraries but allow the user to change if desired
OPTION(BUILD_SHARED_LIBS "Global flag to cause add_library to create shared libraries if on" ON)

# Add to path in order to pick up the FindXXX.cmake files included in this project
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Our code requires a compiler with c++11 support
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

add_definitions(-Wfatal-errors)

# This library depends on ROOT
find_package(ROOT REQUIRED COMPONENTS Table MathMore Geom)

include_directories(
	${ROOT_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}
	include/
)

add_library(tpcrs
	St_base/StMessMgr.cxx
	StarClassLibrary/StParticleTable.cc
	StarClassLibrary/StParticleDefinition.cc
	StarClassLibrary/StPhysicalHelix.cc
	StarClassLibrary/StHelix.cc
	StarClassLibrary/StHyperTriton.cc
	StarClassLibrary/StAntiAlpha.cc
	StarClassLibrary/StAntiHyperTriton.cc
	StarClassLibrary/StAntiDeuteron.cc
	StarClassLibrary/StAntiTriton.cc
	StarClassLibrary/StAntiHelium3.cc
	StarClassLibrary/StHDibaryon.cc
	StarMagField/StarMagField.cxx
	StarRoot/TRArray.cxx
	StarRoot/TRMatrix.cxx
	StarRoot/TRVector.cxx
	StarRoot/TUnixTime.cxx
	StDbUtilities/StGlobalCoordinate.cc
	StDbUtilities/StGlobalDirection.cc
	StDbUtilities/StMagUtilities.cxx
	StDbUtilities/StTpcPadCoordinate.cc
	StDbUtilities/StTpcCoordinate.cxx
	StDbUtilities/StTpcCoordinateTransform.cc
	StTpcRSMaker/StTpcRSMaker.cxx
	StTpcRSMaker/Altro.cxx
	StTpcRSMaker/TF1F.cxx
	StDetectorDbMaker/St_tpcChargeEventC.cxx
	StDetectorDbMaker/St_spaceChargeCorC.cxx
	StEvent/StTpcRawData.cxx
	StEvent/StTpcPixel.cxx
	St_base/StObject.cxx
	StTpcDb/StTpcDb.cxx
)
