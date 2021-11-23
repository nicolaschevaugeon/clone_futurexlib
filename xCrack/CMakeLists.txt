cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xCrack --")
message(STATUS "=============== --")


project(
    xCrack
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

set(XCRACK_SOURCE
# Headers
src/CrackPostpro.h
src/ExportSensors.h
src/FunctionCrackXFEM.h
src/IntegratorCrack.h
src/IntegratorSingularCrack2D.h
src/IntegratorSingularCrack3D.h
src/Jint3D.h
src/JintModalParks.h
src/JintParks.h
src/SingularCrackMapping.h
src/lCrack.h
src/xCrackGeometry.h
src/xGatherMesh.h
src/xcAnalyticalSolutions.h
src/xcCrackBase.h
src/xcFilters.h
src/xcFormLinearEnergyRelease.h
src/xcFront.h
src/xcFrontDomain.h
src/xcFrontDomainManager.h
src/xcFrontManager.h
src/xcFrontPart.h
src/xcFrontSpaceManager.h
src/xcFrontSpaces.h
src/xcInteractionIntegralsOnCrack.h
src/xcValueJs.h
# sources
src/CrackPostpro.cc
src/ExportSensors.cc
src/FunctionCrackXFEM.cc
src/IntegratorCrack.cc
src/IntegratorSingularCrack2D.cc
src/IntegratorSingularCrack3D.cc
src/Jint3D.cc
src/JintModalParks.cc
src/JintParks.cc
src/SingularCrackMapping.cc
src/lCrack.cc
src/xGatherMesh.cc
src/xcAnalyticalSolutions.cc
src/xcFormLinearEnergyRelease.cc
src/xcFront.cc
src/xcFrontDomain.cc
src/xcFrontDomainManager.cc
src/xcFrontManager.cc
src/xcFrontPart.cc
src/xcFrontSpaceManager.cc
src/xcFrontSpaces.cc
src/xcInteractionIntegralsOnCrack.cc
src/xcValueJs.cc
)


add_library(xCrack ${BUILD_SHARED_LIBS} ${XCRACK_SOURCE})
target_include_directories(xCrack PUBLIC src)

target_link_libraries(xCrack PUBLIC  xFEM::xFEM xExport::xExport xCut::xLegacySimpleCut xLinAlg::xInterfaceSuperLu)

# Target properties
set_target_properties(xCrack PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

# Make alias
add_library(xCrack::xCrack ALIAS xCrack)
# Export to a Config.cmake script
export(TARGETS xCrack NAMESPACE xCrack:: FILE xCrackConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xCrack)

# Install
install_target_in_destination(xFEM)

# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
#add_subdirectory(test)
#endif()


# OLD STUFF TO remove !! :
#if(NOT USE_ANN)
#  message(WARNING "Option USE_ANN is OFF but usefull for module xCrack")
#endif()

#if (USE_ANN)
# find_package(ANN REQUIRED)
# add_definitions( -DHAVE_ANN)
# include_directories( ${ANN_INCLUDE_DIR})
#endif(USE_ANN)

#file(GLOB sources src/* )
#add_library(xCrack ${BUILD_SHARED_LIBS}  ${sources}  )
#set_target_properties(xCrack   PROPERTIES COMPILE_FLAGS "-Wno-deprecated ${FLAGS}")


#target_link_libraries(xCrack xLegacySimpleCut )

#install_target_in_destination(xCrack)
