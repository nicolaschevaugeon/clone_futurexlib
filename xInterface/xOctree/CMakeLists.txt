cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xInterface/xOctree --")
message(STATUS "=============== --")


project(
    xInterfacexOctree
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

set(XINTERFACEXOCTREE_SOURCE
# Headers
src/AdaptOctreeToAOMD.h
src/AdaptOctreeToAOMD_imp.h
src/InterfaceOctreeToAOMD.h
src/InterfaceOctreeToAOMD_imp.h
src/OctreeUtilities.h
src/surf2LevelSet.h
# Sources
src/AdaptOctreeToAOMD.cc
src/InterfaceOctreeToAOMD.cc
src/OctreeUtilities.cc
src/surf2LevelSet.cc
)

add_library(xInterfacexOctree ${BUILD_SHARED_LIBS} ${XINTERFACEXOCTREE_SOURCE})
target_include_directories(xInterfacexOctree PUBLIC src)

target_link_libraries(xInterfacexOctree PRIVATE xTool::xTool Trellis::AOMD Trellis::model Trellis::Util xFEM::xFEM xInterface::xInterfaceAOMDGeneral  xLinAlg::xSolverbase PUBLIC xOctree::xOctree)



# Make alias
add_library(xInterface::xInterfacexOctree ALIAS xInterfacexOctree)
# Export to a Config.cmake script
export(TARGETS xInterfacexOctree NAMESPACE xInterface:: FILE xInterfacexOctreeConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfacexOctree)

# Install
install_target_in_destination(xInterfacexOctree)
