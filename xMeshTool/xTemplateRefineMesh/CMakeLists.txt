cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xMeshTools/xTemplateRefineMesh --")
message(STATUS "=============== --")


project(
    xTemplateRefineMesh
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


set(XTEMPLATEREFINEMESH_SOURCE
# Headers   
src/xSplit.h
src/xSplitInfoManager.h
src/xSplitInfoManager_imp.h
src/xSplit_imp.h
)


add_library(xTemplateRefineMesh INTERFACE)
target_include_directories(xTemplateRefineMesh INTERFACE src)

# So that the header files show up in your favorite IDE
target_sources(xTemplateRefineMesh PUBLIC ${XTEMPLATEREFINEMESH_SOURCE})


# Make alias
add_library(xMeshTools::xTemplateRefineMesh ALIAS xTemplateRefineMesh)
# Export to a Config.cmake script
export(TARGETS xTemplateRefineMesh NAMESPACE xMeshTools:: FILE xTemplateRefineMeshConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xTemplateRefineMesh)


# Install
install_target_in_destination(xTemplateRefineMesh)
