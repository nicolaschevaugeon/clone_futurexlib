cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xTensor --")
message(STATUS "=============== --")


project(
    xTensor
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


# Not needed as target xTool has been defined in the same projet !
#find_package(xTool         REQUIRED  )



set(XTENSOR_SOURCE
# Headers   
src/xPoint.h
src/xTensor2.h
src/xTensor2Sym.h
src/xTensor2_imp.h
src/xTensor3.h
src/xTensor4.h
src/xTensor4_imp.h
src/xTensorOperations.h
src/xVector.h
src/xVectorScalarProd.h
src/xVector_imp.h
# Sources
src/xTensor2.cc
)

# Create target
add_library(xTensor ${BUILD_SHARED_LIBS} ${XTENSOR_SOURCE})
# Make alias
add_library(xTensor::xTensor ALIAS xTensor)

target_include_directories(xTensor PUBLIC src)
target_link_libraries(xTensor PUBLIC xTool::xTool)

# Export to a Config.cmake script
export(TARGETS xTensor NAMESPACE xTensor:: FILE xTensorConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xTensor)
