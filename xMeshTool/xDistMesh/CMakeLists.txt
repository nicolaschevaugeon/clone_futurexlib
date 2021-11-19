cmake_minimum_required(VERSION 3.13...3.19)

message(STATUS "=============== --")
message(STATUS "entering xMeshTool/xDistMesh --")
message(STATUS "=============== --")


project(
    xDistMesh
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
#find_package(xTool                 REQUIRED )
#find_package(Trellis                 REQUIRED )
#find_package(xFEM                 REQUIRED )

set(XDISTMESH_SOURCE
# Headers
src/xLoadBalanceTools.h
src/xLoadBalanceTools_imp.h
# Sources
src/xLoadBalanceTools.cc
)

add_library(xDistMesh ${BUILD_SHARED_LIBS} ${XDISTMESH_SOURCE})
target_include_directories(xDistMesh PUBLIC src)
target_link_libraries(xDistMesh PRIVATE xInterface::xInterfaceAOMDGeneral)
target_link_libraries(xDistMesh PUBLIC  xInterface::xInterfaceParMetis)
# Make alias
add_library(xMeshTool::xDistMesh ALIAS xDistMesh)
# Export to a Config.cmake script
export(TARGETS xDistMesh NAMESPACE xMeshTool:: FILE xExportConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xDistMesh)

# Install
install_target_in_destination(xDistMesh)

# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
add_subdirectory(test)
#endif()
