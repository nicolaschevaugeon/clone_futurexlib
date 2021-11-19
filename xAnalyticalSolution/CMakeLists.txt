cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xAnalyticalSolution --")
message(STATUS "=============== --")


project(
    xAnalyticalSolution
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



set(XANALYTICALSOLUTION_SOURCE
# Headers   
Elasticity2D/complexPotential.h
# Sources
Elasticity2D/complexPotential.cc
)


add_library(xAnalyticalSolution ${BUILD_SHARED_LIBS} ${XANALYTICALSOLUTION_SOURCE})
target_include_directories(xAnalyticalSolution PUBLIC Elasticity2D)
target_link_libraries(xAnalyticalSolution PUBLIC xFEM::xFEM)


# Make alias
add_library(xAnalyticalSolution::xAnalyticalSolution ALIAS xAnalyticalSolution)
# Export to a Config.cmake script
export(TARGETS xAnalyticalSolution NAMESPACE xAnalyticalSolution:: FILE xAnalyticalSolutionConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xAnalyticalSolution)

# Install
install_target_in_destination(xAnalyticalSolution)
