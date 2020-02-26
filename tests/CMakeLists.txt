add_executable(test_tpcrs
	test_tpcrs.cxx
	test_tpcrs_dict.cxx
)

root_generate_dictionary(
	test_tpcrs_dict
	tests/GeantEvent.h
	LINKDEF tests/LinkDef.h
)

target_link_libraries(test_tpcrs tpcrs ${ROOT_LIBRARIES})
