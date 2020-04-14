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


function(ADD_TESTS tests)
    foreach(_test_args ${tests})
        # Convert string to list of arguments with the first one being the test name
        string(REGEX MATCHALL "([^\ ]+)" _test_arg_list "${_test_args}")
        list(GET _test_arg_list 0 _name)
        list(GET _test_arg_list 3 _label)
        if($<CONFIG:Debug> OR CMAKE_BUILD_TYPE STREQUAL "Debug")
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
    endforeach()
endfunction()

# name nevents ecutoff label diff_Debug diff_Release
# diff is in lines
set(_tests
"starY16_dAu200    2 0.001  quick 0 0"
"starY15_pp200a    2 0.001  quick 0 0"
"starY15_pp200b    1 0.001  quick 0 0"
"starY14_AuAu200a  1 0.001  quick 0 0"
"starY14_AuAu200b  2 0.0001 quick 0 0"
"starY14_He3Au200  1 0.001  quick 0 0"
"starY12_pp500     1 0.001  quick 0 0"
"starY11_pp500     3 0.0001 quick 0 0"
"starY10_AuAu11    1 0.0001 quick 0 0"
#
"starY16_dAu200   -1 0.001  long  0 0"
"starY15_pp200a   -1 0.001  long  0 0"
"starY15_pp200b   -1 0.001  long  0 1"
"starY14_AuAu200a -1 0.001  long  0 0"
"starY14_AuAu200b -1 0.0001 long  0 0"
"starY14_He3Au200 -1 0.001  long  0 2"
"starY12_CuAu200  -1 0.001  long  0 4"
"starY12_UU200    -1 0.001  long  0 6"
"starY12_pp500    -1 0.001  long  0 0"
"starY11_pp500    -1 0.0001 long  0 0"
"starY10_AuAu62   -1 0.0001 long  0 0"
"starY10_AuAu11   -1 0.0001 long  0 2"
)

add_tests("${_tests}")
