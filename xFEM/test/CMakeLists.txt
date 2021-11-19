cmake_minimum_required(VERSION 3.13...3.19)

message(STATUS "=============== --")
message(STATUS "entering xFEM/test --")
message(STATUS "=============== --")

set(TEST_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/testxOctreeGrid
    ${CMAKE_CURRENT_SOURCE_DIR}/testxRegularGrid
)



create_tests_from_list(LIST_TESTS ${TEST_LIST} TARGETS xFEM::xFEM xInterface::xInterfaceAOMDGeneral)

