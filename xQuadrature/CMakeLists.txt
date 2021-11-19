cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xQuadrature --")
message(STATUS "=============== --")


project(
    xQuadrature
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

# Not needed as these targets have been defined in the same projet !
#find_package(Trellis               REQUIRED )

set(XQUADRATURE_SOURCE
# Headers
src/xGaussLegendreSimplex.h
src/xIntPt.h
src/xIntegrator.h
# Sources
src/xGaussLegendre1D.cc
src/xGaussLegendreSimplex.cc
src/xGaussQuadratureHex.cc
src/xGaussQuadratureQuad.cc
src/xGaussQuadratureTet.cc
src/xGaussQuadratureTri.cc
)

add_library(xQuadrature ${BUILD_SHARED_LIBS} ${XQUADRATURE_SOURCE})
target_include_directories(xQuadrature PUBLIC src)
#target_include_directories(xQuadrature PRIVATE ${Trellis_INCLUDE_DIR})
target_link_libraries(xQuadrature PRIVATE Trellis::AOMD xMapping::xMapping INTERFACE xInterface::xInterfaceAOMDGeneral)


# Make alias
add_library(xQuadrature::xQuadrature ALIAS xQuadrature)
# Export to a Config.cmake script
export(TARGETS xQuadrature NAMESPACE xQuadrature:: FILE xQuadratureConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xQuadrature)

# Install
install_target_in_destination(xQuadrature)
