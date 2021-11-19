cmake_minimum_required(VERSION 3.13...3.19)
message(STATUS "=============== --")
message(STATUS "atomic test: xTool --")
message(STATUS "=============== --")

set(XFILES_CMAKE_UTIL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../xUtil/cmakeUtil"
    CACHE PATH
    "path to Xfiles xUtil/cmakeUtil")

list(APPEND CMAKE_MODULE_PATH "${XFILES_CMAKE_UTIL_PATH}")

include(xutil_ctest_functions)

set(LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/testPartitionManager
    ${CMAKE_CURRENT_SOURCE_DIR}/testDataExchanger
    ${CMAKE_CURRENT_SOURCE_DIR}/testxMPITag
    ${CMAKE_CURRENT_SOURCE_DIR}/testSendOnlyKeysTraits
    ${CMAKE_CURRENT_SOURCE_DIR}/xExportStringDist
    ${CMAKE_CURRENT_SOURCE_DIR}/xDeltaTime
    ${CMAKE_CURRENT_SOURCE_DIR}/xDeltaMemory
     )
create_tests_from_list(LIST_TESTS ${LIST} TARGETS xTool::xTool)

# test with pthread: GREG: NECESSAIRE ???
#list(APPEND EXTERNAL_LIBRARIES   -lpthread )

set(LISTPTHREAD     ${CMAKE_CURRENT_SOURCE_DIR}/mpiInitThreadMultiple      )
create_tests_from_list(LIST_TESTS ${LISTPTHREAD} TARGETS xTool::xTool)
