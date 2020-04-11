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
find_package(ROOT REQUIRED COMPONENTS Cint Table MathMore Geom)

# In case of 32-bit ROOT add the compatible compiler flag
if(${ROOT_CXX_FLAGS} MATCHES "-m32")
	message(STATUS "tpc-rs: Found -m32 option in $ROOT_CXX_FLAGS (root-config). Will proceed with 32 bit build")
	set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS FALSE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
endif()


include(ExternalProject)

# yaml-cpp
ExternalProject_Add(
	yaml-cpp
	URL        https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.3.tar.gz
	URL_HASH   MD5=b45bf1089a382e81f6b661062c10d0c2
	PREFIX     ${CMAKE_BINARY_DIR}/contrib/
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
	           -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=ON
)

set(YAML_CPP_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/contrib/)

add_library(yaml-cpp-lib STATIC IMPORTED)
set_property(TARGET yaml-cpp-lib PROPERTY IMPORTED_LOCATION ${YAML_CPP_INSTALL_PREFIX}/lib/libyaml-cpp.a)
add_dependencies(yaml-cpp-lib yaml-cpp)

set(CMAKE_INSTALL_DATADIR "share"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where read-only data installed")

set(TPCRS_DATA_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/tpc-rs/data"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where project config data is installed")

configure_file(src/config.h.in tpcrs/config.h @ONLY)

add_library(tpcrs
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
	StarClassLibrary/StAlpha.cc
	StarClassLibrary/StAntiBMesonZero.cc
	StarClassLibrary/StAntiBsMesonZero.cc
	StarClassLibrary/StAntiDMesonZero.cc
	StarClassLibrary/StAntiDStarMesonZero.cc
	StarClassLibrary/StAntiKaonZero.cc
	StarClassLibrary/StAntiLambda.cc
	StarClassLibrary/StAntiLambda1520.cc
	StarClassLibrary/StAntiLambdacPlus.cc
	StarClassLibrary/StAntiNeutrinoE.cc
	StarClassLibrary/StAntiNeutrinoMu.cc
	StarClassLibrary/StAntiNeutrinoTau.cc
	StarClassLibrary/StAntiNeutron.cc
	StarClassLibrary/StAntiOmegaMinus.cc
	StarClassLibrary/StAntiOmegacZero.cc
	StarClassLibrary/StAntiProton.cc
	StarClassLibrary/StAntiSigmaMinus.cc
	StarClassLibrary/StAntiSigmaMinus1385.cc
	StarClassLibrary/StAntiSigmaPlus.cc
	StarClassLibrary/StAntiSigmaPlus1385.cc
	StarClassLibrary/StAntiSigmaZero.cc
	StarClassLibrary/StAntiSigmacPlus.cc
	StarClassLibrary/StAntiSigmacPlusPlus.cc
	StarClassLibrary/StAntiSigmacZero.cc
	StarClassLibrary/StAntiXiMinus.cc
	StarClassLibrary/StAntiXiZero.cc
	StarClassLibrary/StAntiXicPlus.cc
	StarClassLibrary/StAntiXicZero.cc
	StarClassLibrary/StBMesonMinus.cc
	StarClassLibrary/StBMesonPlus.cc
	StarClassLibrary/StBMesonZero.cc
	StarClassLibrary/StBsMesonZero.cc
	StarClassLibrary/StCerenkov.cc
	StarClassLibrary/StDMesonMinus.cc
	StarClassLibrary/StDMesonPlus.cc
	StarClassLibrary/StDMesonZero.cc
	StarClassLibrary/StDStarMesonMinus.cc
	StarClassLibrary/StDStarMesonPlus.cc
	StarClassLibrary/StDStarMesonZero.cc
	StarClassLibrary/StDalitz.cc
	StarClassLibrary/StDeuteron.cc
	StarClassLibrary/StDsMesonMinus.cc
	StarClassLibrary/StDsMesonPlus.cc
	StarClassLibrary/StElectron.cc
	StarClassLibrary/StEta.cc
	StarClassLibrary/StEtaPrime.cc
	StarClassLibrary/StGamma.cc
	StarClassLibrary/StGeantino.cc
	StarClassLibrary/StH0Strangelet.cc
	StarClassLibrary/StHe3.cc
	StarClassLibrary/StHelium3.cc
	StarClassLibrary/StJPsi.cc
	StarClassLibrary/StKStarZero.cc
	StarClassLibrary/StKaonMinus.cc
	StarClassLibrary/StKaonPlus.cc
	StarClassLibrary/StKaonZero.cc
	StarClassLibrary/StKaonZeroLong.cc
	StarClassLibrary/StKaonZeroMode0809.cc
	StarClassLibrary/StKaonZeroShort.cc
	StarClassLibrary/StLambda.cc
	StarClassLibrary/StLambda1520.cc
	StarClassLibrary/StLambdacPlus.cc
	StarClassLibrary/StMuonMinus.cc
	StarClassLibrary/StMuonPlus.cc
	StarClassLibrary/StNeutrinoE.cc
	StarClassLibrary/StNeutrinoMu.cc
	StarClassLibrary/StNeutrinoTau.cc
	StarClassLibrary/StNeutron.cc
	StarClassLibrary/StOmegaMeson.cc
	StarClassLibrary/StOmegaMinus.cc
	StarClassLibrary/StOmegacZero.cc
	StarClassLibrary/StOpticalPhoton.cc
	StarClassLibrary/StParticleTable.cc
	StarClassLibrary/StPhi.cc
	StarClassLibrary/StPionMinus.cc
	StarClassLibrary/StPionPlus.cc
	StarClassLibrary/StPionZero.cc
	StarClassLibrary/StPositron.cc
	StarClassLibrary/StProton.cc
	StarClassLibrary/StPsi2s.cc
	StarClassLibrary/StRhoMinus.cc
	StarClassLibrary/StRhoPlus.cc
	StarClassLibrary/StRhoZero.cc
	StarClassLibrary/StSigmaMinus.cc
	StarClassLibrary/StSigmaMinus1385.cc
	StarClassLibrary/StSigmaPlus.cc
	StarClassLibrary/StSigmaPlus1385.cc
	StarClassLibrary/StSigmaZero.cc
	StarClassLibrary/StSigmacPlus.cc
	StarClassLibrary/StSigmacPlusPlus.cc
	StarClassLibrary/StSigmacZero.cc
	StarClassLibrary/StTauMinus.cc
	StarClassLibrary/StTauPlus.cc
	StarClassLibrary/StTriton.cc
	StarClassLibrary/StUpsilon.cc
	StarClassLibrary/StUpsilon2S.cc
	StarClassLibrary/StUpsilon3S.cc
	StarClassLibrary/StWMinusBoson.cc
	StarClassLibrary/StWPlusBoson.cc
	StarClassLibrary/StXiMinus.cc
	StarClassLibrary/StXiZero.cc
	StarClassLibrary/StXiZero1530.cc
	StarClassLibrary/StXicPlus.cc
	StarClassLibrary/StXicZero.cc
	StarClassLibrary/StZZeroBoson.cc
	StarMagField/StarMagField.cxx
	StdEdxY2Maker/StTpcdEdxCorrection.cxx
	StDbUtilities/StGlobalCoordinate.cc
	StDbUtilities/StGlobalDirection.cc
	StDbUtilities/StMagUtilities.cxx
	StDbUtilities/StTpcPadCoordinate.cc
	StDbUtilities/StTpcCoordinate.cxx
	StDbUtilities/StTpcCoordinateTransform.cc
	StTpcRSMaker/StTpcRSMaker.cxx
	StTpcRSMaker/Altro.cxx
	StTpcRSMaker/TF1F.cxx
	StDetectorDbMaker/StDetectorDbChairs.cxx
	StDetectorDbMaker/St_tpcChargeEventC.cxx
	StDetectorDbMaker/St_spaceChargeCorC.cxx
	StTpcDb/StTpcDb.cxx
	src/configurator.cpp
	src/logger.cpp
)

target_include_directories(tpcrs PRIVATE ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}
	include/ ${YAML_CPP_INSTALL_PREFIX}/include/ ${CMAKE_BINARY_DIR})
add_dependencies(tpcrs yaml-cpp-lib)

add_subdirectory(tests)
