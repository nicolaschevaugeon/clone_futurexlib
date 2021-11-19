cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xCut/xLegacySimpleCut --")
message(STATUS "=============== --")


project(
    xLegacySimpleCut
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

set(XLEGACYSIMPLECUT_SOURCE
# Headers
src/xMeshCut.h
src/xPhysSurfByTagging.h
src/xPhysSurfParameter.h
src/xPhysSurfVLS.h
src/xRefCut.h
src/xRefCutToAOMD.h
src/xRefCutToIsoZeroVector.h
src/xRefMesh.h
src/xPhysSurf.h
# Sources
src/xPhysSurfByTagging.cc
src/xPhysSurfParameter.cc
src/xPhysSurfVLS.cc
src/xRefCut.cc
src/xRefCutToAOMD.cc
src/xRefCutToIsoZeroVector.cc
src/xRefMesh.cc
src/xMeshCut.cc
src/xPhysSurf.cc
)

add_library(xLegacySimpleCut ${BUILD_SHARED_LIBS} ${XLEGACYSIMPLECUT_SOURCE})
target_include_directories(xLegacySimpleCut PUBLIC src)


target_link_libraries(xLegacySimpleCut PRIVATE xTool::xTool xFEM::xFEM xInterface::xInterfaceAOMDGeneral)




# Make alias
add_library(xCut::xLegacySimpleCut ALIAS xLegacySimpleCut)
# Export to a Config.cmake script
export(TARGETS xLegacySimpleCut NAMESPACE xCut:: FILE xLegacySimpleCutConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xLegacySimpleCut)


# Install
install_target_in_destination(xLegacySimpleCut)
