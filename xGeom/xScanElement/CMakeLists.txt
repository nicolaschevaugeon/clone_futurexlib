cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xGeom/xScanElement --")
message(STATUS "=============== --")


project(
    xScanElement
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
#find_package(xTensor         REQUIRED  )
#find_package(Trellis         REQUIRED  )

set(XSCANELEMENT_SOURCE
# Headers   
src/xScanElement.h
src/xScanElement_imp.h
# Sources
src/xScanElement.cc
)

add_library(xScanElement ${BUILD_SHARED_LIBS} ${XSCANELEMENT_SOURCE})
target_include_directories(xScanElement PUBLIC src)
target_link_libraries(xScanElement PRIVATE xTensor::xTensor Trellis::AOMD Trellis::model Trellis::Util)


# Make alias
add_library(xGeom::xScanElement ALIAS xScanElement)
# Export to a Config.cmake script
export(TARGETS xScanElement NAMESPACE xGeom::xScanElement FILE xScanElementConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xScanElement)

# Install
install_target_in_destination(xScanElement)
