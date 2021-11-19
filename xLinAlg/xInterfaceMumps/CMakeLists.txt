cmake_minimum_required(VERSION 3.13...3.19)

message(STATUS "=============== --")
message(STATUS "entering xLinAlg/xInterfaceMumps --")
message(STATUS "=============== --")




project(
    xInterfaceMumps
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



find_package(MUMPS REQUIRED COMPONENTS PAR)
find_package(BLAS REQUIRED)

enable_language(Fortran)
find_package(MPI REQUIRED COMPONENTS CXX Fortran)

include(CMakePrintHelpers)
cmake_print_properties(TARGETS MUMPS::MUMPS PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)


#if(SMUMPS_FOUND)
    #message("Interface will use SMUMPS")
##     add_compile_definitions() was introduced in cmake 3.12.
##     add_compile_definitions(USE_SMUMPS)
    #add_definitions(-DUSE_SMUMPS)
#endif()

#if(ZMUMPS_FOUND)
    #message("Interface will use ZMUMPS")
##     add_compile_definitions() was introduced in cmake 3.12.
##     add_compile_definitions(USE_ZMUMPS)
    #add_definitions(-DUSE_ZMUMPS)
#endif()



set(XINTERFACEMUMPS_SOURCE
# Headers
src/xLinearSystemSolverMumps.h
src/xLinearSystemSolverMumpsTraitPolicy.h
src/xLinearSystemSolverMumps_imp.h
# Sources
src/xLinearSystemSolverMumps.cc
)

add_library(xInterfaceMumps   ${BUILD_SHARED_LIBS}  ${XINTERFACEMUMPS_SOURCE}  )
target_include_directories(xInterfaceMumps PUBLIC src)
target_link_libraries(xInterfaceMumps PUBLIC MUMPS::MUMPS xLinAlg::xSolverbase INTERFACE  BLAS::BLAS MPI::MPI_Fortran)
add_dependencies(xInterfaceMumps xLinAlg::xSolverbase) # So that everything is built in the right order !



# Make alias
add_library(xLinAlg::xInterfaceMumps ALIAS xInterfaceMumps)
# Export to a Config.cmake script
export(TARGETS xInterfaceMumps NAMESPACE xLinAlg:: FILE xInterfaceMumpsConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfaceMumps)


# Install
install_target_in_destination(xInterfaceMumps)
