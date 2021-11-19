cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xLinAlg/xInterfaceLapack --")
message(STATUS "=============== --")


project(
    xInterfaceLapack
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


find_package(LAPACK        REQUIRED)
# Not needed as these targets have been defined in the same projet !
#find_package(xSolverBase   REQUIRED)

set(XINTERFACELAPACK_SOURCE
# Headers
src/xLapackInterface.h
# Sources
src/xLapackInterface.cc
)

add_library(xInterfaceLapack ${BUILD_SHARED_LIBS} ${XINTERFACELAPACK_SOURCE})
target_include_directories(xInterfaceLapack PUBLIC src)
#target_include_directories(xInterfaceLapack PRIVATE ${xSolverBase_INCLUDE_DIR})
target_link_libraries(xInterfaceLapack PUBLIC xLinAlg::xSolverbase PRIVATE LAPACK::LAPACK)


# Make alias
add_library(xLinAlg::xInterfaceLapack ALIAS xInterfaceLapack)
# Export to a Config.cmake script
export(TARGETS xInterfaceLapack NAMESPACE xLinAlg:: FILE xInterfaceLapackConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfaceLapack)


# Install
install_target_in_destination(xInterfaceLapack)

# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
add_subdirectory(test)
#endif()
