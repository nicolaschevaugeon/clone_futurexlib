cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xInterface/xInterfaceParMetis --")
message(STATUS "=============== --")


project(
    xInterfaceParMetis
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

find_package(METIS REQUIRED)

message(STATUS "PARMETIS_VERSION :${PARMETIS_VERSION}")
message(STATUS "PARMETIS_MAJOR_VERSION :${PARMETIS_MAJOR_VERSION}")


set(XINTERFACEPARMETIS_SOURCE
# Headers
src/ParMetisInterface.h
src/ParMetisInterface_types.h
# Sources
src/ParMetisInterface.cc
)

add_library(xInterfaceParMetis ${BUILD_SHARED_LIBS} ${XINTERFACEPARMETIS_SOURCE})
target_include_directories(xInterfaceParMetis PUBLIC src)

target_link_libraries(xInterfaceParMetis PRIVATE xInterface::xInterfaceAOMDGeneral)

target_link_libraries(xInterfaceParMetis PUBLIC  METIS::METIS)


# Make alias
add_library(xInterface::xInterfaceParMetis ALIAS xInterfaceParMetis)
# Export to a Config.cmake script
export(TARGETS xInterfaceParMetis NAMESPACE xInterface:: FILE xExportConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xInterfaceParMetis)

# Install
install_target_in_destination(xInterfaceParMetis)
