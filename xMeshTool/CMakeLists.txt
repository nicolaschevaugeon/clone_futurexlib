cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xMeshTools --")
message(STATUS "=============== --")


project(
    xMeshTools
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



#if(BUILD_XDISTMESH)
    #add_subdirectory(xDistMesh)
#endif(BUILD_XDISTMESH)

add_subdirectory(xReader)

