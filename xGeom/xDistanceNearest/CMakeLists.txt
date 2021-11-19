cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xDistanceNearest --")
message(STATUS "=============== --")


project(
    xDistanceNearest
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


set(XDISTANCENEAREST_SOURCE
# Headers   
src/xDistanceNearestPoint.h
src/xDistanceNearestPointGenerator.h
src/xDistanceNearestPointGenerator_imp.h
src/xDistanceNearestPoint_internal.h
src/xDistanceNearestPoint_internal_imp.h
src/xDistanceNearestPoint_traits.h
# Sources
src/xDistanceNearestPointGenerator.cc
)

# Create target
add_library(xDistanceNearest ${BUILD_SHARED_LIBS} ${XDISTANCENEAREST_SOURCE})
target_include_directories(xDistanceNearest PUBLIC src)
target_link_libraries(xDistanceNearest PRIVATE Trellis::AOMD Trellis::model Trellis::Util xInterface::xInterfaceAOMDGeneral xFEM::xFEM)


# Make alias
add_library(xGeom::xDistanceNearest ALIAS xDistanceNearest)
# Export to a Config.cmake script
export(TARGETS xDistanceNearest NAMESPACE xGeom:: FILE xGeomToolsConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xDistanceNearest)


# Install
install_target_in_destination(xDistanceNearest)


