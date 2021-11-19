cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xInterface/xTemplateRefineMesh --")
message(STATUS "=============== --")


project(
    xInterfacexTemplateRefineMesh
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

set(XTEMPLATEREFINEMESH_SOURCE
# Headers
src/xSplitAOMD.h
src/xSplitAOMD_imp.h
# Sources
src/xSplitAOMD.cc
)

add_library(xInterfacexTemplateRefineMesh ${BUILD_SHARED_LIBS} ${XTEMPLATEREFINEMESH_SOURCE})
target_include_directories(xInterfacexTemplateRefineMesh PUBLIC src)


target_link_libraries(xInterfacexTemplateRefineMesh PRIVATE xFEM::xFEM xInterface::xInterfaceAOMDGeneral PUBLIC xMeshTools::xTemplateRefineMesh)




# Make alias
add_library(xInterface::xTemplateRefineMesh ALIAS xInterfacexTemplateRefineMesh)
# Export to a Config.cmake script
export(TARGETS xInterfacexTemplateRefineMesh NAMESPACE xInterface:: FILE xInterfacexTemplateRefineMeshConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfacexTemplateRefineMesh)

# Install
install_target_in_destination(xInterfacexTemplateRefineMesh)
