cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xFEM   --")
message(STATUS "=============== --")


project(
    xFEM
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



set(XFEM_SOURCE
# Headers
src/xExchange.h
src/xAlgorithm.h
src/xApplyCommandOnIntegrationRule.h
src/xApproxFctPtr.h
src/xApproxFunction.h
src/xApproxFunctionBySupportComponent.h
src/xApproxFunctionCorrectedShifted.h
src/xApproxFunctionDerivDiscXFEMOctree.h
src/xApproxFunctionExtendedShape.h
src/xApproxFunctionHighOrder.h
src/xApproxFunctionHighOrderQH.h
src/xAssembler.h
src/xAssemblerDispatcher.h
src/xAssemblerDispatcherTraitPolicy.h
src/xAssemblerDispatcher_imp.h
src/xAssemblerTraitPolicy.h
src/xCommandOnGeomElem.h
src/xData.h
src/xDataExchangerTools.h
src/xDebug.h
src/xDispatcher.h
src/xDistStateOfValue.h
src/xDistStateOfValue_imp.h
src/xDomain.h
src/xElement.h
src/xElementAccess.h
src/xElementAccess_imp.h
src/xEntityFilter.h
src/xEntityToEntity.h
src/xEnv.h
src/xEval.h
src/xEvalAspectRatio.h
src/xEvalStorage.h
src/xExchange_imp.h
src/xSubMesh.h
src/xExtendShapeFcts.h
src/xExtendShapeFcts_imp.h
src/xFactory.h
src/xFemMatrix.h
src/xField.h
src/xFieldStorage.h
src/xField_imp.h
src/xFiniteElement.h
src/xForm.h
src/xFormOptimized.h
src/xGeomElem.h
src/xGetSupport.h
src/xIntegrationRule.h
src/xIntegrationRuleStored.h
src/xIter.h
src/xLagrangeMappingBuilder.h
src/xLevelSet.h
src/xLevelSetOperators.h
src/xLinkOnFrontLinkGenerator.h
src/xMappingBuilder.h
src/xMappingBuilderHolder.h
src/xMaterial.h
src/xMaterialManager.h
src/xMaterialSensitivity.h
src/xMesh.h
src/xNearestNeighborInterface.h
src/xNonLocalInfoForKeysAndFcts.h
src/xNonLocalInfoForKeysAndFcts_imp.h
src/xOperators.h
src/xParseData.h
src/xPartition.h
src/xPointToDouble.h
src/xRegion.h
src/xRenumbering.h
src/xRenumbering_imp.h
src/xSolVisitor.h
src/xSpace.h
src/xSpaceExtendedShape.h
src/xSpaceFactory.h
src/xSpacePtr.h
src/xSpaceFactory_imp.h
src/xSpaceMultiXFEM.h
src/xSpacePolynomial.h
src/xSpacePolynomialQH.h
src/xStateOfValue.h
src/xStateOfValue_imp.h
src/xSubMeshManager.h
src/xSupportComponent.h
src/xTable.h
src/xTag.h
src/xTensors.h
src/xTensorsPtr.h
src/xValKey.h
src/xValKeyDataManager.h
src/xValManager.h
src/xValue.h
src/xValueCreatorLinkOnFront.h
src/xValueCreatorLinkOnFrontFiltered.h
src/xValueCreatorLinkOnFrontFiltered_imp.h
src/xValueCreatorLinkOnFront_imp.h
src/xValueCreators.h
src/xValueLinearCombination.h
src/xValueLinearCombination_imp.h
src/xValueManager.h
src/xValue_imp.h
src/xVariabManager.h
src/xVectorField.h
src/xVectorLevelSet.h
src/xVectorLevelSet_imp.h
src/xZone.h
# Sources
src/xApproxFunction.cc
src/xApproxFunctionBySupportComponent.cc
src/xApproxFunctionCorrectedShifted.cc
src/xApproxFunctionExtendedShape.cc
src/xApproxFunctionHighOrder.cc
src/xApproxFunctionHighOrderBernsteinQH.cc
src/xApproxFunctionHighOrderHierarchicalLegendreQH.cc
src/xApproxFunctionHighOrderLagrangeQH.cc
src/xApproxFunctionHighOrderQH.cc
src/xCommandOnGeomElem.cc
src/xRegion.cc
src/xDataExchangerTools.cc
src/xDebug.cc
src/xElement.cc
src/xEntityFilter.cc
src/xEntityToEntity.cc
src/xEnv.cc
src/xEvalStorage.cc
src/xFemMatrix.cc
src/xFiniteElement.cc
src/xForm.cc
src/xGeomElem.cc
src/xGetSupport.cc
src/xIntegrationRule.cc
src/xIntegrationRuleStored.cc
src/xLagrangeMappingBuilder.cc
src/xLevelSet.cc
src/xLevelSetOperators.cc
src/xLinkOnFrontLinkGenerator.cc
src/xMappingBuilderHolder.cc
src/xData.cc
src/xMaterial.cc
src/xMaterialManager.cc
src/xMesh.cc
src/xNearestNeighborInterface.cc
src/xNonLocalInfoForKeysAndFcts.cc
src/xAlgorithm.cc
src/xParseData.cc
src/xSubMesh.cc
src/xSpace.cc
src/xSpaceExtendedShape.cc
src/xSpaceFactory.cc
src/xSpaceMultiXFEM.cc
src/xSpacePolynomial.cc
src/xSpacePolynomialQH.cc
src/xStateOfValue.cc
src/xSubMeshManager.cc
src/xSupportComponent.cc
src/xTensors.cc
src/xValKey.cc
src/xValue.cc
src/xValueCreatorLinkOnFrontFiltered.cc
src/xValueCreators.cc
src/xValueLinearCombination.cc
src/xVectorField.cc
src/xVectorLevelSet.cc
src/xZone.cc
)


add_library(xFEM ${BUILD_SHARED_LIBS} ${XFEM_SOURCE})
target_include_directories(xFEM PUBLIC src)



# GREG: BBBEAAAARRRRKKK : a corriger ! (mais un peu mieux maintenant : INTERFACE_INCLUDE_DIRECTORIES ne donne que la partie publique des includes)
#get_target_property(XINTERFACEAOMDGENERAL_INCLUDES xInterface::xInterfaceAOMDGeneral INTERFACE_INCLUDE_DIRECTORIES)
#get_target_property(XSOLVERBASE_INCLUDES xLinAlg::xSolverbase INTERFACE_INCLUDE_DIRECTORIES)
#get_target_property(XINTERFACELAPACK_INCLUDES xLinAlg::xInterfaceLapack INTERFACE_INCLUDE_DIRECTORIES)
#target_include_directories(xFEM PRIVATE ${XINTERFACEAOMDGENERAL_INCLUDES} ${XSOLVERBASE_INCLUDES} ${XINTERFACELAPACK_INCLUDES})

target_link_libraries(xFEM PUBLIC xTool::xTool Trellis::AOMD Trellis::model Trellis::Util xTensor::xTensor xMapping::xMapping xGeom::xGeomTools xQuadrature::xQuadrature xMeshTools::xReader xLinAlg::xSolverbase xLinAlg::xInterfaceLapack)


# Target properties
set_target_properties(xFEM PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")


# Make alias
add_library(xFEM::xFEM ALIAS xFEM)
# Export to a Config.cmake script
export(TARGETS xFEM NAMESPACE xFEM:: FILE xFEMConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xFEM)

# Install
install_target_in_destination(xFEM)

# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
add_subdirectory(test)
#endif()
