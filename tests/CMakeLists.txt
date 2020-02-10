add_executable(test_tpcrs
	test_tpcrs.cxx
)

target_link_libraries(test_tpcrs
	tpcrs ${ROOT_LIBRARIES}
)
