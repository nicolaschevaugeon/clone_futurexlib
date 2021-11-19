cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xLinAlg --")
message(STATUS "=============== --")


project(
    xLinAlg
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



#---------------------------------------------------------
#----------------------- SOLVERBASE ----------------------
#---------------------------------------------------------
add_subdirectory(xSolverBase)


#---------------------------------------------------------
#----------------------- INTERFACES ----------------------
#---------------------------------------------------------

#if(BUILD_XINTERFACE_SUPERLU)
  add_subdirectory(xInterfaceSuperLu)
#endif(BUILD_XINTERFACE_SUPERLU)


#if(BUILD_XINTERFACE_LAPACK)
  add_subdirectory(xInterfaceLapack)
#endif(BUILD_XINTERFACE_LAPACK)


#if(BUILD_XINTERFACE_MUMPS)
  add_subdirectory(xInterfaceMumps)
#endif(BUILD_XINTERFACE_MUMPS)

