cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering Trellis --")
message(STATUS "=============== --")


project(
    Trellis
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

#########################################
# Model library
#########################################
add_library(model ${BUILD_SHARED_LIBS} "")
add_subdirectory(model)


#########################################
# Util library
#########################################
add_library(Util ${BUILD_SHARED_LIBS} "")
add_subdirectory(Util)


#########################################
# AOMD library
#########################################
add_library(AOMD ${BUILD_SHARED_LIBS} "")
add_subdirectory(AOMD)
# AOMD depends on model and Util targets (libraries)
target_link_libraries(AOMD PRIVATE model Util)

# Make alias
add_library(Trellis::model ALIAS model)
add_library(Trellis::Util ALIAS Util)
add_library(Trellis::AOMD ALIAS AOMD)
# Export to a Config.cmake script
export(TARGETS model Util AOMD NAMESPACE xGeom FILE TrellisConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE Trellis)


# Install
install_target_in_destination(model Util AOMD)
