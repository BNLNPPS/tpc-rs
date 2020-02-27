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
	StDetectorDbMaker/StDetectorDbChairs.cxx
	StDetectorDbMaker/St_tpcChargeEventC.cxx
	StDetectorDbMaker/St_spaceChargeCorC.cxx
	StEvent/StTpcRawData.cxx
	StEvent/StTpcPixel.cxx
	St_base/StObject.cxx
	St_db_Maker/St_db_Maker.cxx
	St_db_Maker/StValiSet.cxx
	StTpcDb/StTpcDb.cxx
	StChain/StMaker.cxx
	src/tpcrs_tables.cxx
	tpcrs_tables_dict.cxx
)

root_generate_dictionary(
	tpcrs_tables_dict
	tables/St_asic_thresholds_Table.h
	tables/St_beamInfo_Table.h
	tables/St_g2t_tpc_hit_Table.h
	tables/St_g2t_track_Table.h
	tables/St_g2t_vertex_Table.h
	tables/St_itpcPadGainT0_Table.h
	tables/St_iTPCSurvey_Table.h
	tables/St_MagFactor_Table.h
	tables/St_MDFCorrection_Table.h
	tables/St_richvoltages_Table.h
	tables/St_spaceChargeCor_Table.h
	tables/St_starClockOnl_Table.h
	tables/St_starMagOnl_Table.h
	tables/St_Survey_Table.h
	tables/St_tpcAltroParams_Table.h
	tables/St_tpcAnodeHV_Table.h
	tables/St_tpcAnodeHVavg_Table.h
	tables/St_TpcAvgCurrent_Table.h
	tables/St_TpcAvgPowerSupply_Table.h
	tables/St_TpcEffectivedX_Table.h
	tables/St_TpcSecRowCor_Table.h
	tables/St_tpcCalibResolutions_Table.h
	tables/St_tpcChargeEvent_Table.h
	tables/St_tpcCorrection_Table.h
	tables/St_tpcDimensions_Table.h
	tables/St_tpcDriftVelocity_Table.h
	tables/St_tpcEffectiveGeom_Table.h
	tables/St_tpcElectronics_Table.h
	tables/St_tpcFieldCage_Table.h
	tables/St_tpcFieldCageShort_Table.h
	tables/St_tpcGas_Table.h
	tables/St_tpcGlobalPosition_Table.h
	tables/St_tpcGridLeak_Table.h
	tables/St_tpcHighVoltages_Table.h
	tables/St_tpcHVPlanes_Table.h
	tables/St_tpcOmegaTau_Table.h
	tables/St_tpcPadConfig_Table.h
	tables/St_tpcPadGainT0_Table.h
	tables/St_tpcPadPlanes_Table.h
	tables/St_tpcPadResponse_Table.h
	tables/St_tpcPadrowT0_Table.h
	tables/St_tpcPedestal_Table.h
	tables/St_tpcRDOMap_Table.h
	tables/St_tpcRDOMasks_Table.h
	tables/St_tpcRDOT0offset_Table.h
	tables/St_TpcResponseSimulator_Table.h
	tables/St_tpcSCGL_Table.h
	tables/St_tpcSectorPosition_Table.h
	tables/St_tpcSectorT0offset_Table.h
	tables/St_tpcSlowControlSim_Table.h
	tables/St_tpcWirePlanes_Table.h
	tables/St_trgTimeOffset_Table.h
	tables/St_trigDetSums_Table.h
	tables/St_tss_tsspar_Table.h
	StDetectorDbMaker/St_beamInfoC.h
	StDetectorDbMaker/St_MagFactorC.h
	StDetectorDbMaker/St_MDFCorrectionC.h
	StDetectorDbMaker/St_TpcAdcCorrectionBC.h
	StDetectorDbMaker/St_TpcAdcCorrectionMDF.h
	StDetectorDbMaker/St_TpcCurrentCorrectionC.h
	StDetectorDbMaker/St_TpcdChargeC.h
	StDetectorDbMaker/St_TpcdEdxCorC.h
	StDetectorDbMaker/St_TpcDriftDistOxygenC.h
	StDetectorDbMaker/St_TpcdXCorrectionBC.h
	StDetectorDbMaker/St_TpcEdgeC.h
	StDetectorDbMaker/St_TpcEffectivedXC.h
	StDetectorDbMaker/St_tpcGasTemperatureC.h
	StDetectorDbMaker/St_TpcLengthCorrectionBC.h
	StDetectorDbMaker/St_TpcLengthCorrectionMDF.h
	StDetectorDbMaker/St_tpcMethaneInC.h
	StDetectorDbMaker/St_TpcMultiplicityC.h
	StDetectorDbMaker/St_TpcPadCorrectionMDF.h
	StDetectorDbMaker/St_TpcPhiDirectionC.h
	StDetectorDbMaker/St_tpcPressureBC.h
	StDetectorDbMaker/St_TpcrChargeC.h
	StDetectorDbMaker/St_TpcRowQC.h
	StDetectorDbMaker/St_TpcSecRowBC.h
	StDetectorDbMaker/St_TpcSecRowCC.h
	StDetectorDbMaker/St_TpcSecRowCorC.h
	StDetectorDbMaker/St_TpcSpaceChargeC.h
	StDetectorDbMaker/St_TpcTanLC.h
	StDetectorDbMaker/St_tpcTimeDependenceC.h
	StDetectorDbMaker/St_tpcWaterOutC.h
	StDetectorDbMaker/St_TpcZCorrectionBC.h
	StDetectorDbMaker/St_TpcZDCC.h
	LINKDEF include/LinkDef.h
)

add_subdirectory(tests)
