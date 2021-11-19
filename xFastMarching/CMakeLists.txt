cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xFastMarching --")
message(STATUS "=============== --")


project(
    xFastMarching
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



set(XFASTMARCHING_SOURCE
# Headers   
src/FM.h
src/FMDataExchanger.h
src/FMDist.h
src/FMDistInternal.h
src/FMEntityStorage.h
src/FMInternal.h
src/FMUpdater.h
src/FMUpdaterDist.h
src/FMUtil.h
src/Skeleton/FMSK_fastmarching.h
src/Skeleton/FMSK_fastmarching_updater.h
src/Skeleton/FMSK_meshutil.h
src/Skeleton/FMSK_skeleton_enriched_levelset.h
src/Skeleton/FM_seq.h
src/linearalgebra3d.h
src/linearalgebra3dPackUnPack.h
src/meshinterfaceAOMD.h
src/meshinterfacexHM.h
src/meshinterfacexRegion.h
)

# Create target
add_library(xFastMarching INTERFACE)
target_include_directories(xFastMarching INTERFACE src)

# So that the header files show up in your favorite IDE
target_sources(xFastMarching PUBLIC ${XFASTMARCHING_SOURCE})

# Make alias
add_library(xFastMarching::xFastMarching ALIAS xFastMarching)

# Export to a Config.cmake script
export(TARGETS xFastMarching NAMESPACE xFastMarching:: FILE xFastMarchingConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xFastMarching)

# Install
install_target_in_destination(xFastMarching)

