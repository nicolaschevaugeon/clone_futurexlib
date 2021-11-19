cmake_minimum_required(VERSION 3.13...3.19)
message(STATUS "=============== --")
message(STATUS "atomic test: xMapping --")
message(STATUS "=============== --")

enable_testing()

set(TEST_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/testxLagrangeMapping
     )

create_tests_from_list(LIST_TESTS ${TEST_LIST} TARGETS xMapping::xMapping xGeom::xGeomTools)



