cmake_minimum_required(VERSION 3.10)
project(xMeshInterface)
message(STATUS "---------------------------")
message(STATUS "PROJECT: xMeshInterface  --")
message(STATUS "---------------------------")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../Util/cmakeUtil/  )
include(common_functions)
set(BUILD_SHARED_LIBS "SHARED")
define_archos_suffixe(ARCHOS)


find_package(xTool)
find_package(xTensor)
find_package(xGeomTools)
find_package(xFEM)
 

message(STATUS "PROJECT: xMeshInterface Base --")
file(GLOB sources    ${CMAKE_CURRENT_SOURCE_DIR}/miBase/*.cc  )

set( INCLUDES_DIR  
  ${CMAKE_CURRENT_SOURCE_DIR}/miBase
  ${MPI_INCLUDE_PATH}
  ${xTool_INCLUDE_DIR}
  ${xTensor_INCLUDE_DIR}
  ${xGeomTools_INCLUDE_DIR}
  ${xFEM_INCLUDE_DIR} 
)
include_directories(${INCLUDES_DIR} )
add_library(miBase ${BUILD_SHARED_LIBS}  ${sources}  )
set_target_properties(miBase   PROPERTIES COMPILE_FLAGS "-Wno-deprecated ${FLAGS}"  )
install(TARGETS miBase DESTINATION ${DEVROOT}/lib/${ARCHOS})



find_package(Trellis)
find_package(xInterfaceAOMDGeneral)
find_package(xTemplateRefineMesh)



message(STATUS "PROJECT: xMeshInterface AOMD --")
find_package(Trellis)
file(GLOB sources    ${CMAKE_CURRENT_SOURCE_DIR}/miAOMD/*.cc  )
set( INCLUDES_DIR  
  ${CMAKE_CURRENT_SOURCE_DIR}/miAOMD
  ${CMAKE_CURRENT_SOURCE_DIR}/miBase
  ${xTools_INCLUDE_DIR}
  ${xFEM_INCLUDE_DIR} 
  ${xTensor_INCLUDE_DIR}
  ${Trellis_INCLUDE_DIR} 
  ${xTemplateRefineMesh_INCLUDE_DIR}
  ${xInterfaceAOMDGeneral_INCLUDE_DIR}
  ${xFEM_INCLUDE_DIR} 
  ${MPI_INCLUDE_PATH}
  )


include_directories(${INCLUDES_DIR} )
add_library(miAOMD ${BUILD_SHARED_LIBS}  ${sources}  )
set_target_properties(miAOMD   PROPERTIES COMPILE_FLAGS "-Wno-deprecated ${FLAGS}"  )
install(TARGETS miAOMD DESTINATION ${DEVROOT}/lib/${ARCHOS})
