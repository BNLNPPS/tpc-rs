cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project(tpc-rs LANGUAGES CXX)

set(TPCRS_VERSION_MAJOR 0 CACHE STRING "major version of tpc-rs" FORCE)
set(TPCRS_VERSION_MINOR 0 CACHE STRING "minor version of tpc-rs" FORCE)
set(TPCRS_VERSION_PATCH 7 CACHE STRING "patch version of tpc-rs" FORCE)
set(TPCRS_VERSION "${TPCRS_VERSION_MAJOR}.${TPCRS_VERSION_MINOR}.${TPCRS_VERSION_PATCH}" CACHE STRING "version of tpc-rs" FORCE)

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

# Set a debug postfix
set(CMAKE_DEBUG_POSTFIX "-dbg")

# Set installation destinations
if(UNIX)
    include(GNUInstallDirs)

    set(TPCRS_INC_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/tpcrs")
    set(TPCRS_RUNTIME_INSTALL_DIR "${CMAKE_INSTALL_BINDIR}")
    set(TPCRS_LIBRARY_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
    set(TPCRS_ARCHIVE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
    set(TPCRS_FRAMEWORK_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")

    set(TPCRS_ADDITIONAL_FILES_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/tpc-rs-${TPCRS_VERSION}")
    set(TPCRS_CMAKE_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/tpc-rs-${TPCRS_VERSION}/cmake")
    set(TPCRS_DATA_INSTALL_DIR "${CMAKE_INSTALL_FULL_DATADIR}/tpc-rs-${TPCRS_VERSION}/data")
elseif(WIN32)
    set(TPCRS_INC_INSTALL_DIR       "include/tpcrs")
    set(TPCRS_RUNTIME_INSTALL_DIR   "bin")
    set(TPCRS_LIBRARY_INSTALL_DIR   "bin")
    set(TPCRS_ARCHIVE_INSTALL_DIR   "lib")
    set(TPCRS_FRAMEWORK_INSTALL_DIR "bin")

    set(TPCRS_ADDITIONAL_FILES_INSTALL_DIR "share/tpc-rs-${TPCRS_VERSION}")
    set(TPCRS_CMAKE_CONFIG_INSTALL_DIR "share/tpc-rs-${TPCRS_VERSION}/cmake")
    set(TPCRS_DATA_INSTALL_DIR "share/tpc-rs-${TPCRS_VERSION}/data")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    message(FATAL_ERROR "Could not set install folders for this platform!")
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

add_subdirectory(src)
add_subdirectory(tests)

# Install additional files such as README, LICENSE, etc.
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md" DESTINATION ${TPCRS_ADDITIONAL_FILES_INSTALL_DIR})

install(EXPORT tpcrsTargets DESTINATION ${TPCRS_CMAKE_CONFIG_INSTALL_DIR} FILE tpcrs-config.cmake)
