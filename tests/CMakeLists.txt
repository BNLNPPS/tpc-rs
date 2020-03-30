add_executable(test_tpcrs
	test_tpcrs.cxx
	test_tpcrs_dict.cxx
)

root_generate_dictionary(
	test_tpcrs_dict
	tests/GeantEvent.h
	LINKDEF tests/LinkDef.h
	OPTIONS -I${CMAKE_SOURCE_DIR}/include/
	INCDIRS ${CMAKE_SOURCE_DIR}
)

target_include_directories(test_tpcrs PRIVATE ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include/ ${YAML_CPP_INSTALL_PREFIX}/include/)
target_link_libraries(test_tpcrs tpcrs ${ROOT_LIBRARIES} ${YAML_CPP_INSTALL_PREFIX}/lib/libyaml-cpp.a)
