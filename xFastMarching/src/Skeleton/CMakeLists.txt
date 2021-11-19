message(" FMSkeleton ")

find_package(xTool REQUIRED)
find_package(xFEM REQUIRED)
find_package(xFastMarching REQUIRED)
find_package(xExport REQUIRED)
find_package(Trellis REQUIRED)
find_package(xInterfaceAOMDGeneral REQUIRED)
find_package(xSolverBase REQUIRED)
find_package(xInterfaceLapack REQUIRED)
find_package(BLAS2 REQUIRED)


list(APPEND EXTERNAL_INCLUDES
  ${MPI_CXX_INCLUDE_PATH} 
  ${xTool_INCLUDE_DIR} 
  ${Trellis_INCLUDE_DIR}  
  ${xFEM_INCLUDE_DIR}
  ${xFastMarching_INCLUDE_DIR}
  ${xExport_INCLUDE_DIR}
  ${xInterfaceAOMDGeneral_INCLUDE_DIR}
  )

message("##### ${BLAS_LIBRARIES}")
message("##### ${LAPACK_LIBRARIES}")


list(APPEND EXTERNAL_LIBRARIES 
  ${xFEM_LIBRARIES}
  ${xExport_LIBRARIES}
  ${xSolverBase_LIBRARIES}
  ${xInterfaceLapack_LIBRARIES}
  ${Lapack_LIBRARIES}
  ${Trellis_LIBRARIES}
  ${xInterfaceAOMDGeneral_LIBRARIES}
  ${xTool_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
 # "blas"
#  "/usr/lib/liblapack.so.3"
  )


add_executable("test_skeleton" main_sk.cc data)
include_directories(                ${EXTERNAL_INCLUDES} )
target_link_libraries("test_skeleton" ${EXTERNAL_LIBRARIES})

set_target_properties("test_skeleton"
  PROPERTIES  
  COMPILE_FLAGS "-std=c++11 -Wno-deprecated -Wall -fmax-errors=3"
  )


add_custom_command(OUTPUT data COMMAND ln -s  ${CMAKE_CURRENT_SOURCE_DIR}/data)

