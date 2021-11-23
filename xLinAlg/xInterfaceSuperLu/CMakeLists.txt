cmake_minimum_required(VERSION 3.13...3.19)

message(STATUS "=============== --")
message(STATUS "entering xLinAlg/xInterfaceSuperLu --")
message(STATUS "=============== --")




project(
    xInterfaceSuperLu
    VERSION 1.0
    LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Enforce C++17 as required by the compiler
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    
    # By default, we build shared libs
    set(BUILD_SHARED_LIBS "SHARED")
endif()




find_package(SuperLU REQUIRED)
string(REGEX MATCH "[0-9]+" SuperLU_MAJOR_VERSION "${SuperLU_VERSION}")

configure_file (
    src/SuperLU_Version.h.in
    "${CMAKE_CURRENT_BINARY_DIR}/SuperLU_Version.h"
)

find_package(BLAS REQUIRED)


include(CMakePrintHelpers)
cmake_print_properties(TARGETS SuperLU::SuperLU PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)


set(XINTERFACESUPERLU_SOURCE
# Headers
"${CMAKE_CURRENT_BINARY_DIR}/SuperLU_Version.h"
src/xLinearSystemSolverSuperLU.h
src/xLinearSystemSolverSuperLUDataType.h
src/xLinearSystemSolverSuperLUTraitPolicy.h
src/xLinearSystemSolverSuperLU_imp.h
src/xLinearSystemSolverSuperLUcplx.h
# Sources
src/xLinearSystemSolverSuperLU.cc
)

add_library(xInterfaceSuperLu   ${BUILD_SHARED_LIBS}  ${XINTERFACESUPERLU_SOURCE}  )
target_include_directories(xInterfaceSuperLu PUBLIC src "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(xInterfaceSuperLu PUBLIC SuperLU::SuperLU xLinAlg::xSolverBase PRIVATE BLAS::BLAS)
add_dependencies(xInterfaceSuperLu xLinAlg::xSolverBase) # So that everything is built in the right order !



# Make alias
add_library(xLinAlg::xInterfaceSuperLu ALIAS xInterfaceSuperLu)
# Export to a Config.cmake script
export(TARGETS xInterfaceSuperLu NAMESPACE xLinAlg:: FILE xInterfaceSuperLuConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfaceSuperLu)


# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
#add_subdirectory(test)
#endif()
