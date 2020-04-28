cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project(TPCRS LANGUAGES CXX)

include(CTest)

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

# This library depends on ROOT
find_package(ROOT REQUIRED COMPONENTS Cint Table MathMore Geom)

# In case of 32-bit ROOT add the compatible compiler flag
if(${ROOT_CXX_FLAGS} MATCHES "-m32")
    message(STATUS "tpc-rs: Found -m32 option in $ROOT_CXX_FLAGS (root-config). Will proceed with 32 bit build")
    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS FALSE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
    set(CMAKE_SIZEOF_VOID_P 4)
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

# Setup yaml-cpp library as external project
include(cmake/ExternalProject_yaml-cpp.cmake)

add_subdirectory(src)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Install additional files such as README, LICENSE, etc.
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md" DESTINATION ${TPCRS_ADDITIONAL_FILES_INSTALL_DIR})

install(EXPORT tpcrsTargets DESTINATION ${TPCRS_CMAKE_CONFIG_INSTALL_DIR} FILE tpcrs-config.cmake)
