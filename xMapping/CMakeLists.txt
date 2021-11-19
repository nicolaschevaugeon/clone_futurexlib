cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xMapping --")
message(STATUS "=============== --")


project(
    xMapping
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
#find_package(xTool         REQUIRED)
#find_package(xTensor       REQUIRED)
#find_package(xGeomTools    REQUIRED)


set(XMAPPING_SOURCE
# Headers  
src/xLagrangeMapping.h
src/xMapping.h
src/xReferenceElement.h
# Sources
src/xLagrangeMapping.cc
src/xMapping.cc
)


add_library(xMapping ${BUILD_SHARED_LIBS} ${XMAPPING_SOURCE})
target_include_directories(xMapping PUBLIC src)
target_link_libraries(xMapping PRIVATE xTool::xTool xTensor::xTensor xGeom::xGeomTools)
add_dependencies(xMapping xTensor xTool xGeomTools) # So that everything is built in the right order
#target_link_libraries(xMapping PRIVATE xTool::xTool xTensor::xTensor xGeom::xGeomTools xInterface::xInterfaceAOMDGeneral)
#add_dependencies(xMapping xTensor xTool xGeomTools xInterfaceAOMDGeneral) # So that everything is built in the right order

# Make alias
add_library(xMapping::xMapping ALIAS xMapping)
# Export to a Config.cmake script
export(TARGETS xMapping NAMESPACE xMapping:: FILE xMappingConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xMapping)

# Install
install_target_in_destination(xMapping)



# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
add_subdirectory(test)
#endif()
