cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xExport --")
message(STATUS "=============== --")


project(
    xExport
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


set(XEXPORT_SOURCE
# Headers
src/xExport.h
src/xExportAlgorithm.h
src/xExportEnsight.h
src/xExportGmsh.h
src/xExportGmsh_imp.h
src/xExportManager.h
src/xExportTranslate.h
src/xPostPro.h
src/xPostProMSH.h
src/xPostProVTU.h
src/xPostProcessingManager.h
src/xSensors.h
# Sources
src/xExport.cc
src/xExportEnsight.cc
src/xExportManager.cc
src/xExportTranslate.cc
src/xPostPro.cc
src/xPostProMSH.cc
src/xPostProVTU.cc
src/xPostProcessingManager.cc
src/xSensors.cc
src/xExportAlgorithm.cc
src/xExportGmsh.cc
)

add_library(xExport ${BUILD_SHARED_LIBS} ${XEXPORT_SOURCE})
target_include_directories(xExport PUBLIC src)


target_link_libraries(xExport PUBLIC xFEM::xFEM PRIVATE xInterface::xInterfaceAOMDGeneral)




# Make alias
add_library(xExport::xExport ALIAS xExport)
# Export to a Config.cmake script
export(TARGETS xExport NAMESPACE xExport:: FILE xExportConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xExport)


# Install
install_target_in_destination(xExport)
