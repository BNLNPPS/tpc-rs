configure_file(config.h.in config.h @ONLY)

set(header_files
    ${TPCRS_SOURCE_DIR}/include/tpcrs/configurator.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/config_type.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/config_yaml.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/tpcrs.h
    ${TPCRS_SOURCE_DIR}/include/tpcrs/detail/simulator.h
)

set(source_files
    particles/StParticleTable.cc
    particles/StParticleDefinition.cc
    particles/StHyperTriton.cc
    particles/StAntiAlpha.cc
    particles/StAntiHyperTriton.cc
    particles/StAntiDeuteron.cc
    particles/StAntiTriton.cc
    particles/StAntiHelium3.cc
    particles/StHDibaryon.cc
    particles/StAlpha.cc
    particles/StAntiBMesonZero.cc
    particles/StAntiBsMesonZero.cc
    particles/StAntiDMesonZero.cc
    particles/StAntiDStarMesonZero.cc
    particles/StAntiKaonZero.cc
    particles/StAntiLambda.cc
    particles/StAntiLambda1520.cc
    particles/StAntiLambdacPlus.cc
    particles/StAntiNeutrinoE.cc
    particles/StAntiNeutrinoMu.cc
    particles/StAntiNeutrinoTau.cc
    particles/StAntiNeutron.cc
    particles/StAntiOmegaMinus.cc
    particles/StAntiOmegacZero.cc
    particles/StAntiProton.cc
    particles/StAntiSigmaMinus.cc
    particles/StAntiSigmaMinus1385.cc
    particles/StAntiSigmaPlus.cc
    particles/StAntiSigmaPlus1385.cc
    particles/StAntiSigmaZero.cc
    particles/StAntiSigmacPlus.cc
    particles/StAntiSigmacPlusPlus.cc
    particles/StAntiSigmacZero.cc
    particles/StAntiXiMinus.cc
    particles/StAntiXiZero.cc
    particles/StAntiXicPlus.cc
    particles/StAntiXicZero.cc
    particles/StBMesonMinus.cc
    particles/StBMesonPlus.cc
    particles/StBMesonZero.cc
    particles/StBsMesonZero.cc
    particles/StCerenkov.cc
    particles/StDMesonMinus.cc
    particles/StDMesonPlus.cc
    particles/StDMesonZero.cc
    particles/StDStarMesonMinus.cc
    particles/StDStarMesonPlus.cc
    particles/StDStarMesonZero.cc
    particles/StDalitz.cc
    particles/StDeuteron.cc
    particles/StDsMesonMinus.cc
    particles/StDsMesonPlus.cc
    particles/StElectron.cc
    particles/StEta.cc
    particles/StEtaPrime.cc
    particles/StGamma.cc
    particles/StGeantino.cc
    particles/StH0Strangelet.cc
    particles/StHe3.cc
    particles/StHelium3.cc
    particles/StJPsi.cc
    particles/StKStarZero.cc
    particles/StKaonMinus.cc
    particles/StKaonPlus.cc
    particles/StKaonZero.cc
    particles/StKaonZeroLong.cc
    particles/StKaonZeroMode0809.cc
    particles/StKaonZeroShort.cc
    particles/StLambda.cc
    particles/StLambda1520.cc
    particles/StLambdacPlus.cc
    particles/StMuonMinus.cc
    particles/StMuonPlus.cc
    particles/StNeutrinoE.cc
    particles/StNeutrinoMu.cc
    particles/StNeutrinoTau.cc
    particles/StNeutron.cc
    particles/StOmegaMeson.cc
    particles/StOmegaMinus.cc
    particles/StOmegacZero.cc
    particles/StOpticalPhoton.cc
    particles/StParticleTable.cc
    particles/StPhi.cc
    particles/StPionMinus.cc
    particles/StPionPlus.cc
    particles/StPionZero.cc
    particles/StPositron.cc
    particles/StProton.cc
    particles/StPsi2s.cc
    particles/StRhoMinus.cc
    particles/StRhoPlus.cc
    particles/StRhoZero.cc
    particles/StSigmaMinus.cc
    particles/StSigmaMinus1385.cc
    particles/StSigmaPlus.cc
    particles/StSigmaPlus1385.cc
    particles/StSigmaZero.cc
    particles/StSigmacPlus.cc
    particles/StSigmacPlusPlus.cc
    particles/StSigmacZero.cc
    particles/StTauMinus.cc
    particles/StTauPlus.cc
    particles/StTriton.cc
    particles/StUpsilon.cc
    particles/StUpsilon2S.cc
    particles/StUpsilon3S.cc
    particles/StWMinusBoson.cc
    particles/StWPlusBoson.cc
    particles/StXiMinus.cc
    particles/StXiZero.cc
    particles/StXiZero1530.cc
    particles/StXicPlus.cc
    particles/StXicZero.cc
    particles/StZZeroBoson.cc
    coords.cpp
    mag_field.cpp
    mag_field_utils.cpp
    dedx_correction.cpp
    dedx_parameterization.cpp
    simulator.cpp
    altro.cpp
    TF1F.cpp
    struct_containers.cpp
    configurator.cpp
    logger.cpp
    math_cephes.cpp
    math_funcs.cpp
    track_helix.cpp
)

add_library(tpcrs ${header_files} ${source_files})

target_include_directories(tpcrs
    INTERFACE
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:include/tpcrs/contrib>
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${YAML_CPP_INSTALL_PREFIX}/include ${CMAKE_BINARY_DIR}/src ${ROOT_INCLUDE_DIR})

target_link_libraries(tpcrs PRIVATE yaml-cpp-lib ${ROOT_LIBRARIES})

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
    ${TPCRS_SOURCE_DIR}/data/P10T.root
    ${TPCRS_SOURCE_DIR}/data/bfield_full_negative_2D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_full_negative_3D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_full_positive_2D.dat
    ${TPCRS_SOURCE_DIR}/data/bfield_full_positive_3D.dat
    ${TPCRS_SOURCE_DIR}/data/const_full_positive_2D.dat
    ${TPCRS_SOURCE_DIR}/data/dNdE_Bichsel.root
    ${TPCRS_SOURCE_DIR}/data/dNdx_Bichsel.root
    ${TPCRS_SOURCE_DIR}/data/dNdx_Heed.root
)

install(DIRECTORY ${YAML_CPP_INSTALL_PREFIX}/include/yaml-cpp DESTINATION ${TPCRS_INC_INSTALL_DIR}/contrib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION ${TPCRS_INC_INSTALL_DIR})
install(FILES ${header_files} DESTINATION ${TPCRS_INC_INSTALL_DIR})
install(FILES ${data_files}   DESTINATION ${TPCRS_DATA_INSTALL_DIR})
install(FILES ${version_file} DESTINATION ${TPCRS_CMAKE_CONFIG_INSTALL_DIR})
