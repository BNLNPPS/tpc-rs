add_executable(test_tpcrs
    test_tpcrs.cpp
    test_tpcrs_dict.cxx
)

root_generate_dictionary(
    test_tpcrs_dict
    include/tpcrs/structs.h
    tests/GeantEvent.h
    LINKDEF tests/LinkDef.h
    OPTIONS -I${CMAKE_SOURCE_DIR}/include/
    INCDIRS ${CMAKE_SOURCE_DIR}
)

target_include_directories(test_tpcrs PRIVATE ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include/ ${YAML_CPP_INSTALL_PREFIX}/include/)
target_link_libraries(test_tpcrs tpcrs ${ROOT_LIBRARIES} ${YAML_CPP_INSTALL_PREFIX}/lib/libyaml-cpp.a)


include(ExternalProject)

if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(TEST_DATA_FILE 1CJASQerc5KqCyOSEU043VkOfJx_OtI2Z)
else()
    set(TEST_DATA_FILE 1XESdyqg6kXhu5gPrgS3drq3UIBu0f1l5)
endif()

ExternalProject_Add(
    test-data
    DOWNLOAD_COMMAND ${CMAKE_SOURCE_DIR}/utils/gdrive.py ${TEST_DATA_FILE} > test-data.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/contrib
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -E tar xz <DOWNLOAD_DIR>/test-data.tar.gz
    INSTALL_COMMAND ${CMAKE_COMMAND} -E rename <BINARY_DIR>/data ${CMAKE_BINARY_DIR}/tests/data
)


function(ADD_TESTS tests)
    foreach(_test_args ${tests})
        # Convert string to list of arguments with the first one being the test name
        string(REGEX MATCHALL "([^\ ]+)" _test_arg_list "${_test_args}")
        list(GET _test_arg_list 0 _name)
        list(GET _test_arg_list 3 _label)

        if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
            list(GET _test_arg_list 6 _diff_lines)
        elseif(${CMAKE_SIZEOF_VOID_P} EQUAL 4 AND ($<CONFIG:Debug> OR CMAKE_BUILD_TYPE STREQUAL "Debug"))
            list(GET _test_arg_list 4 _diff_lines)
        else()
            list(GET _test_arg_list 5 _diff_lines)
        endif()

        set(_test_name "${_name}_${_label}")
        set(_cmd "command time -f \"time: real %E user %U sys %S\" ./test_tpcrs ${_test_args} &> log_${_test_name}")
        set(_cmd "${_cmd} && diff -y --suppress-common-lines ${_name}_inp.log ${_name}_out.log")
        set(_cmd "${_cmd} | wc -l | (read N; ((\"$N\" == ${_diff_lines})))")

        add_test(NAME "${_test_name}" COMMAND bash -c "${_cmd}")
        set_tests_properties("${_test_name}" PROPERTIES LABELS "${_label}" TIMEOUT 3600)
        set_tests_properties("${_test_name}" PROPERTIES DEPENDS test-data)
    endforeach()
endfunction()

# name nevents ecutoff label diff-32-Debug diff-32-Release diff-64
# diff is in lines
set(_tests
"starY16_dAu200    2 0.001  quick 0 0  0"
"starY15_pp200a    2 0.001  quick 0 0  0"
"starY15_pp200b    1 0.001  quick 0 0  0"
"starY14_AuAu200a  1 0.001  quick 0 0  0"
"starY14_AuAu200b  2 0.0001 quick 0 0  0"
"starY14_He3Au200  1 0.001  quick 0 0  0"
"starY12_pp500     1 0.001  quick 0 0  0"
"starY11_pp500     3 0.0001 quick 0 0  0"
"starY10_AuAu11    1 0.0001 quick 0 0  0"
#
"starY16_dAu200   -1 0.001  long  0 0  0"
"starY15_pp200a   -1 0.001  long  0 0  0"
"starY15_pp200b   -1 0.001  long  0 1 16"
"starY14_AuAu200a -1 0.001  long  0 0 17"
"starY14_AuAu200b -1 0.0001 long  0 0  6"
"starY14_He3Au200 -1 0.001  long  0 2  3"
"starY12_CuAu200  -1 0.001  long  0 4 29"
"starY12_UU200    -1 0.001  long  0 6 24"
"starY12_pp500    -1 0.001  long  0 0  3"
"starY11_pp500    -1 0.0001 long  0 0  0"
"starY10_AuAu62   -1 0.0001 long  0 0 17"
"starY10_AuAu11   -1 0.0001 long  0 2  7"
)

add_tests("${_tests}")
