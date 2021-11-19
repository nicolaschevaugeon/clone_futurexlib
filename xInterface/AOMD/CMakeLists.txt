cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xInterface/AOMD --")
message(STATUS "=============== --")


project(
    xInterfaceAOMDGeneral
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

set(XINTERFACEAOMDGENERAL_SOURCE
# Headers
general/src/xAOMDEntityUtil.h
general/src/xAspectRatioAOMD.h
general/src/xAttachableChar.h
general/src/xAttachableGP.h
general/src/xAttachedDataManagerAOMD.h
general/src/xMeshToolAOMD.h
general/src/xMshToAOMDReader.h
# Sources
#general/src/xAOMDEntityUtil.cc
general/src/xAspectRatioAOMD.cc
general/src/xMeshToolAOMD.cc
general/src/xMshToAOMDReader.cc
)

add_library(xInterfaceAOMDGeneral ${BUILD_SHARED_LIBS} ${XINTERFACEAOMDGENERAL_SOURCE})
target_include_directories(xInterfaceAOMDGeneral PUBLIC general/src)

# GREG: BBBEAAAARRRRKKK : Add xFEM/src as PRIVATE include dir (a corriger !)
target_include_directories(xInterfaceAOMDGeneral PRIVATE $<TARGET_PROPERTY:xFEM,INCLUDE_DIRECTORIES>)


# To display generator result (at build time)
#add_custom_command(TARGET xInterfaceAOMDGeneral POST_BUILD
  #COMMAND ${CMAKE_COMMAND} -E echo 
  #"target dir = $<TARGET_PROPERTY:xFEM,INCLUDE_DIRECTORIES>")


target_link_libraries(xInterfaceAOMDGeneral PUBLIC xTool::xTool Trellis::AOMD Trellis::model Trellis::Util xTensor::xTensor xMapping::xMapping xGeom::xGeomTools xMeshTools::xReader)
target_link_libraries(xInterfaceAOMDGeneral PUBLIC xTool::xTool Trellis::AOMD Trellis::model Trellis::Util xTensor::xTensor xMapping::xMapping xGeom::xGeomTools xMeshTools::xReader)



# Make alias
add_library(xInterface::xInterfaceAOMDGeneral ALIAS xInterfaceAOMDGeneral)
# Export to a Config.cmake script
export(TARGETS xInterfaceAOMDGeneral NAMESPACE xInterface:: FILE xInterfaceAOMDGeneralConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfaceAOMDGeneral)

# Install
install_target_in_destination(xInterfaceAOMDGeneral)
