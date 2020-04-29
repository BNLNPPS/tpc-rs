add_executable(test_tpcrs
    test_tpcrs.cpp
    test_tpcrs_dict.cxx
)

root_generate_dictionary(
    test_tpcrs_dict
    include/g2t_tpc_hit.h
    include/g2t_track.h
    include/g2t_vertex.h
    tests/GeantEvent.h
    LINKDEF tests/LinkDef.h
    OPTIONS -I${CMAKE_SOURCE_DIR}/include/
    INCDIRS ${CMAKE_SOURCE_DIR}
)

target_include_directories(test_tpcrs PRIVATE ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include/ ${YAML_CPP_INSTALL_PREFIX}/include/)
target_link_libraries(test_tpcrs tpcrs ${ROOT_LIBRARIES} ${YAML_CPP_INSTALL_PREFIX}/lib/libyaml-cpp.a)

set(_tests
"starY16_dAu200   2       "
"starY15_pp200a   2       "
"starY15_pp200b   1       "
"starY14_AuAu200a 1       "
"starY14_AuAu200b 2 0.0001"
"starY14_He3Au200 1       "
"starY12_CuAu200  1       "
"starY12_UU200    1       "
"starY12_pp500    1       "
"starY11_pp500    3 0.0001"
"starY10_AuAu62   1 0.0001"
"starY10_AuAu11   1 0.0001"
)

set(_disabled_tests
"starY12_CuAu200"
"starY12_UU200"
"starY10_AuAu62"
)

foreach(_test_args ${_tests})
    # Convert string to list of arguments with the first one being the test name
    string(REGEX MATCHALL "([^\ ]+)" _test_arg_list "${_test_args}")
    list(GET _test_arg_list 0 _name)
    set(_cmd "command time -f \"time: real %E user %U sys %S\" ./test_tpcrs ${_test_args} &> log_${_name} && diff ${_name}_inp.log ${_name}_out.log &> /dev/null")
    set(_disabled FALSE)
    if (${_name} IN_LIST _disabled_tests)
       set(_disabled TRUE)
       continue()
    endif()
    add_test(NAME ${_name} COMMAND bash -c "${_cmd}")
    set_tests_properties(${_name} PROPERTIES TIMEOUT 600 WILL_FAIL ${_disabled})
endforeach()
