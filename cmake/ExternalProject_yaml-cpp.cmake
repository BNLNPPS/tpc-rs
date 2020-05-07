include(ExternalProject)

set(_YAML_CPP_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# yaml-cpp
ExternalProject_Add(
    yaml-cpp
    URL        https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.3.tar.gz
    URL_HASH   MD5=b45bf1089a382e81f6b661062c10d0c2
    PREFIX     ${CMAKE_BINARY_DIR}/contrib/
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_CXX_FLAGS=${_YAML_CPP_CXX_FLAGS}
               -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
)

set(YAML_CPP_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/contrib/)

add_library(yaml-cpp-lib STATIC IMPORTED)
set_property(TARGET yaml-cpp-lib PROPERTY IMPORTED_LOCATION ${YAML_CPP_INSTALL_PREFIX}/lib/libyaml-cpp.a)
add_dependencies(yaml-cpp-lib yaml-cpp)
