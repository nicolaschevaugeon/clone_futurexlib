cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xTool --")
message(STATUS "=============== --")


project(
    xTool
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


find_package(MPI REQUIRED)

set(XTOOL_SOURCE
# Headers   
src/workInProgress.h
src/xAspectRatio.h
src/xDataExchanger.h
src/xDataExchangerTraits.h
src/xDataExchanger_imp.h
src/xDataType.h
src/xDataType_imp.h
src/xDeltaMemory.h
src/xDeltaTime.h
src/xExportStringDist.h
src/xGeneralUnorderedMapDataManager.h
src/xGenericOperations.h
src/xIteratorTools.h
src/xMPIDataType.h
src/xMPIEnv.h
src/xMPITag.h
src/xMacro.h
src/xMemoryMonitor.h
src/xPartitionManager.h
src/xPartitionManagerTools.h
src/xPartitionManagerTools_imp.h
src/xRawDataExchanger.h
src/xSingleton.h
src/xStringManager.h
src/xUnorderedMapDataManager.h
# Sources
src/workInProgress.cc
src/xDeltaMemory.cc
src/xDeltaTime.cc
src/xExportStringDist.cc
src/xMPIDataType.cc
src/xMPIEnv.cc
src/xMPITag.cc
)

# Create target
add_library(xTool ${BUILD_SHARED_LIBS} ${XTOOL_SOURCE})
# Make alias
add_library(xTool::xTool ALIAS xTool)

target_include_directories(xTool PUBLIC src)
target_link_libraries(xTool PUBLIC MPI::MPI_CXX)

# Export to a Config.cmake script
export(TARGETS xTool NAMESPACE xTool:: FILE xToolConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xTool)

# Tests
set(BUILD_ATOMIC_TEST ON
    CACHE BOOL
    "build atomic test for the lib")

if(BUILD_ATOMIC_TEST)
  enable_testing()
  add_subdirectory(test)
endif()

