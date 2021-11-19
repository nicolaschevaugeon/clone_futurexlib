cmake_minimum_required(VERSION 3.13...3.19)
message(STATUS "=============== --")
message(STATUS "atomic test: DISTMESH --")
message(STATUS "=============== --")

 
enable_testing()

set(LISTAOMD   ${CMAKE_CURRENT_SOURCE_DIR}/testPartitionManagerAOMD       )



create_tests_from_list(LIST_TESTS ${LISTAOMD} TARGETS xMeshTool::xDistMesh xTool::xTool xInterface::xInterfaceAOMDGeneral xFEM::xFEM)


#find_package(METIS       REQUIRED )
set(LISTAOMDPARMETIS  ${CMAKE_CURRENT_SOURCE_DIR}/testLoadBalanceToolsAOMD  )
create_tests_from_list(LIST_TESTS ${LISTAOMDPARMETIS} TARGETS xMeshTool::xDistMesh xTool::xTool xInterface::xInterfaceAOMDGeneral xInterface::xInterfaceParMetis xFEM::xFEM)
