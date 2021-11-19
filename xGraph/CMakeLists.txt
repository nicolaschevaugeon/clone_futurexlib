cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xGraph --")
message(STATUS "=============== --")


project(
    xGraph
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



set(XGRAPH_SOURCE
# Headers   
weightedLimitedGraph/src/nodeAndConnectedEdge.h
weightedLimitedGraph/src/nodeAndConnectedEdgeDist.h
weightedLimitedGraph/src/nodeAndConnectedEdgeDistInternal.h
)

# Create target
add_library(xGraph INTERFACE)
target_include_directories(xGraph INTERFACE weightedLimitedGraph/src)

# So that the header files show up in your favorite IDE
target_sources(xGraph PUBLIC ${XGRAPH_SOURCE})


# Make alias
add_library(xGraph::xGraph ALIAS xGraph)
# Export to a Config.cmake script
export(TARGETS xGraph NAMESPACE xGraph:: FILE xGraphConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xGraph)


# Install
install_target_in_destination(xGraph)
