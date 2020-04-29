configure_file(config.h.in config.h @ONLY)

set(header_files
    ${TPCRS_SOURCE_DIR}/include/tpcrs/configurator.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/config_structs.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/config_yaml.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/digi_data.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/enums.h
)

set(source_files
    ../StarClassLibrary/StParticleTable.cc
    ../StarClassLibrary/StParticleDefinition.cc
    ../StarClassLibrary/StPhysicalHelix.cc
    ../StarClassLibrary/StHelix.cc
    ../StarClassLibrary/StHyperTriton.cc
    ../StarClassLibrary/StAntiAlpha.cc
    ../StarClassLibrary/StAntiHyperTriton.cc
    ../StarClassLibrary/StAntiDeuteron.cc
    ../StarClassLibrary/StAntiTriton.cc
    ../StarClassLibrary/StAntiHelium3.cc
    ../StarClassLibrary/StHDibaryon.cc
    ../StarClassLibrary/StAlpha.cc
    ../StarClassLibrary/StAntiBMesonZero.cc
    ../StarClassLibrary/StAntiBsMesonZero.cc
    ../StarClassLibrary/StAntiDMesonZero.cc
    ../StarClassLibrary/StAntiDStarMesonZero.cc
    ../StarClassLibrary/StAntiKaonZero.cc
    ../StarClassLibrary/StAntiLambda.cc
    ../StarClassLibrary/StAntiLambda1520.cc
    ../StarClassLibrary/StAntiLambdacPlus.cc
    ../StarClassLibrary/StAntiNeutrinoE.cc
    ../StarClassLibrary/StAntiNeutrinoMu.cc
    ../StarClassLibrary/StAntiNeutrinoTau.cc
    ../StarClassLibrary/StAntiNeutron.cc
    ../StarClassLibrary/StAntiOmegaMinus.cc
    ../StarClassLibrary/StAntiOmegacZero.cc
    ../StarClassLibrary/StAntiProton.cc
    ../StarClassLibrary/StAntiSigmaMinus.cc
    ../StarClassLibrary/StAntiSigmaMinus1385.cc
    ../StarClassLibrary/StAntiSigmaPlus.cc
    ../StarClassLibrary/StAntiSigmaPlus1385.cc
    ../StarClassLibrary/StAntiSigmaZero.cc
    ../StarClassLibrary/StAntiSigmacPlus.cc
    ../StarClassLibrary/StAntiSigmacPlusPlus.cc
    ../StarClassLibrary/StAntiSigmacZero.cc
    ../StarClassLibrary/StAntiXiMinus.cc
    ../StarClassLibrary/StAntiXiZero.cc
    ../StarClassLibrary/StAntiXicPlus.cc
    ../StarClassLibrary/StAntiXicZero.cc
    ../StarClassLibrary/StBMesonMinus.cc
    ../StarClassLibrary/StBMesonPlus.cc
    ../StarClassLibrary/StBMesonZero.cc
    ../StarClassLibrary/StBsMesonZero.cc
    ../StarClassLibrary/StCerenkov.cc
    ../StarClassLibrary/StDMesonMinus.cc
    ../StarClassLibrary/StDMesonPlus.cc
    ../StarClassLibrary/StDMesonZero.cc
    ../StarClassLibrary/StDStarMesonMinus.cc
    ../StarClassLibrary/StDStarMesonPlus.cc
    ../StarClassLibrary/StDStarMesonZero.cc
    ../StarClassLibrary/StDalitz.cc
    ../StarClassLibrary/StDeuteron.cc
    ../StarClassLibrary/StDsMesonMinus.cc
    ../StarClassLibrary/StDsMesonPlus.cc
    ../StarClassLibrary/StElectron.cc
    ../StarClassLibrary/StEta.cc
    ../StarClassLibrary/StEtaPrime.cc
    ../StarClassLibrary/StGamma.cc
    ../StarClassLibrary/StGeantino.cc
    ../StarClassLibrary/StH0Strangelet.cc
    ../StarClassLibrary/StHe3.cc
    ../StarClassLibrary/StHelium3.cc
    ../StarClassLibrary/StJPsi.cc
    ../StarClassLibrary/StKStarZero.cc
    ../StarClassLibrary/StKaonMinus.cc
    ../StarClassLibrary/StKaonPlus.cc
    ../StarClassLibrary/StKaonZero.cc
    ../StarClassLibrary/StKaonZeroLong.cc
    ../StarClassLibrary/StKaonZeroMode0809.cc
    ../StarClassLibrary/StKaonZeroShort.cc
    ../StarClassLibrary/StLambda.cc
    ../StarClassLibrary/StLambda1520.cc
    ../StarClassLibrary/StLambdacPlus.cc
    ../StarClassLibrary/StMuonMinus.cc
    ../StarClassLibrary/StMuonPlus.cc
    ../StarClassLibrary/StNeutrinoE.cc
    ../StarClassLibrary/StNeutrinoMu.cc
    ../StarClassLibrary/StNeutrinoTau.cc
    ../StarClassLibrary/StNeutron.cc
    ../StarClassLibrary/StOmegaMeson.cc
    ../StarClassLibrary/StOmegaMinus.cc
    ../StarClassLibrary/StOmegacZero.cc
    ../StarClassLibrary/StOpticalPhoton.cc
    ../StarClassLibrary/StParticleTable.cc
    ../StarClassLibrary/StPhi.cc
    ../StarClassLibrary/StPionMinus.cc
    ../StarClassLibrary/StPionPlus.cc
    ../StarClassLibrary/StPionZero.cc
    ../StarClassLibrary/StPositron.cc
    ../StarClassLibrary/StProton.cc
    ../StarClassLibrary/StPsi2s.cc
    ../StarClassLibrary/StRhoMinus.cc
    ../StarClassLibrary/StRhoPlus.cc
    ../StarClassLibrary/StRhoZero.cc
    ../StarClassLibrary/StSigmaMinus.cc
    ../StarClassLibrary/StSigmaMinus1385.cc
    ../StarClassLibrary/StSigmaPlus.cc
    ../StarClassLibrary/StSigmaPlus1385.cc
    ../StarClassLibrary/StSigmaZero.cc
    ../StarClassLibrary/StSigmacPlus.cc
    ../StarClassLibrary/StSigmacPlusPlus.cc
    ../StarClassLibrary/StSigmacZero.cc
    ../StarClassLibrary/StTauMinus.cc
    ../StarClassLibrary/StTauPlus.cc
    ../StarClassLibrary/StTriton.cc
    ../StarClassLibrary/StUpsilon.cc
    ../StarClassLibrary/StUpsilon2S.cc
    ../StarClassLibrary/StUpsilon3S.cc
    ../StarClassLibrary/StWMinusBoson.cc
    ../StarClassLibrary/StWPlusBoson.cc
    ../StarClassLibrary/StXiMinus.cc
    ../StarClassLibrary/StXiZero.cc
    ../StarClassLibrary/StXiZero1530.cc
    ../StarClassLibrary/StXicPlus.cc
    ../StarClassLibrary/StXicZero.cc
    ../StarClassLibrary/StZZeroBoson.cc
    ../StarMagField/StarMagField.cxx
    ../StdEdxY2Maker/StTpcdEdxCorrection.cxx
    ../StDbUtilities/StGlobalCoordinate.cc
    ../StDbUtilities/StGlobalDirection.cc
    ../StDbUtilities/StMagUtilities.cxx
    ../StDbUtilities/StTpcPadCoordinate.cc
    ../StDbUtilities/StTpcCoordinate.cxx
    ../StDbUtilities/StTpcCoordinateTransform.cc
    ../StTpcRSMaker/StTpcRSMaker.cxx
    ../StTpcRSMaker/Altro.cxx
    ../StTpcRSMaker/TF1F.cxx
    ../StDetectorDbMaker/StDetectorDbChairs.cxx
    ../StDetectorDbMaker/St_tpcChargeEventC.cxx
    ../StDetectorDbMaker/St_spaceChargeCorC.cxx
    ../StTpcDb/StTpcDb.cxx
    configurator.cpp
    logger.cpp
    math_cephes.cpp
    math_funcs.cpp
)

add_library(tpcrs ${header_files} ${source_files})

target_include_directories(tpcrs PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include ${YAML_CPP_INSTALL_PREFIX}/include ${CMAKE_BINARY_DIR}/src ${ROOT_INCLUDE_DIR})
add_dependencies(tpcrs yaml-cpp-lib)

set_target_properties(tpcrs PROPERTIES
    VERSION ${TPCRS_VERSION}
    POSITION_INDEPENDENT_CODE ON)

install(TARGETS tpcrs EXPORT tpcrsTargets
    RUNTIME   DESTINATION ${TPCRS_RUNTIME_INSTALL_DIR}
    LIBRARY   DESTINATION ${TPCRS_LIBRARY_INSTALL_DIR}
    ARCHIVE   DESTINATION ${TPCRS_ARCHIVE_INSTALL_DIR}
    FRAMEWORK DESTINATION ${TPCRS_FRAMEWORK_INSTALL_DIR})

# Create and install version file
include(CMakePackageConfigHelpers)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/cmake/tpcrs-config-version.cmake")
write_basic_package_version_file(${version_file}
    VERSION ${TPCRS_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(data_files
    ${TPCRS_SOURCE_DIR}/data/BichselT.root
    ${TPCRS_SOURCE_DIR}/data/P10T.root
    ${TPCRS_SOURCE_DIR}/data/PaiT.root
    ${TPCRS_SOURCE_DIR}/data/StarFieldZ.root
    ${TPCRS_SOURCE_DIR}/data/bfield_full_negative_2D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_full_negative_3D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_full_positive_2D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_full_positive_3D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_half_positive_2D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_half_positive_3D.dat
    ${TPCRS_SOURCE_DIR}/data/const_full_positive_2D.dat
    ${TPCRS_SOURCE_DIR}/data/dNdE_Bichsel.root
    ${TPCRS_SOURCE_DIR}/data/dNdx_Bichsel.root
    ${TPCRS_SOURCE_DIR}/data/dNdx_Heed.root
    ${TPCRS_SOURCE_DIR}/data/steel_magfieldmap.dat
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION ${TPCRS_INC_INSTALL_DIR})
install(FILES ${header_files} DESTINATION ${TPCRS_INC_INSTALL_DIR})
install(FILES ${data_files}   DESTINATION ${TPCRS_DATA_INSTALL_DIR})
install(FILES ${version_file} DESTINATION ${TPCRS_CMAKE_CONFIG_INSTALL_DIR})
