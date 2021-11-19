cmake_minimum_required(VERSION 2.6)
project(xMeshInterfaceAtomicTest)
message(STATUS "----------------------------------------")
message(STATUS "PROJECT: xMeshInterface Atomic Tests  --")
message(STATUS "----------------------------------------")

list(APPEND CMAKE_MODULE_PATH ${DEVROOT}/Util/cmakeUtil/  )
include(common_functions)
set(BUILD_SHARED_LIBS "SHARED")
define_archos_suffixe(ARCHOS)


find_package(xTool           )
find_package(xTensor         )
find_package(Trellis         )
find_package(miBase          )
find_package(miAOMD          )

 
set(LISTTEST
	${CMAKE_CURRENT_SOURCE_DIR}/attachableData
	)

set(EXTERNAL_INCLUDES
	${xTool_INCLUDE_DIR}
	${miBase_INCLUDE_DIR}
	${miAOMD_INCLUDE_DIR}
	${Trellis_INCLUDE_DIR}
	${MPI_CXX_INCLUDE_PATH}
	${xTensor_INCLUDE_DIR}
	)

set(EXTERNAL_LIBRARIES
	${xTool_LIBRARIES}
	${Trellis_LIBRARIES}
	${miBase_LIBRARIES}
	${miAOMD_LIBRARIES}
	${xTensor_LIBRARIES}
	${MPI_CXX_LIBRARIES}
	${xTensor_LIBRARIES}
	)
#create_tests_from_list(${LISTTEST})
create_tests()
