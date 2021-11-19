cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xGeomTools --")
message(STATUS "=============== --")


project(
    xGeomTools
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


set(XGEOMTOOLS_SOURCE
# Headers   
src/xOctreeGrid.h
src/xRegularGrid.h
src/xSimpleGeometry.h
# Sources
src/xBoundingBox.cc
src/xRegularGrid.cc
)

# Create target
add_library(xGeomTools ${BUILD_SHARED_LIBS} ${XGEOMTOOLS_SOURCE})
target_include_directories(xGeomTools PUBLIC src)
target_link_libraries(xGeomTools PRIVATE xTool::xTool PUBLIC xTensor::xTensor)


# Make alias
add_library(xGeom::xGeomTools ALIAS xGeomTools)
# Export to a Config.cmake script
export(TARGETS xGeomTools NAMESPACE xGeom:: FILE xGeomToolsConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xGeomTools)


# Install
install_target_in_destination(xGeomTools)
