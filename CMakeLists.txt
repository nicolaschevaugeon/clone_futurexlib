message(STATUS "=============== --")
message(STATUS "entering Xfiles --")
message(STATUS "=============== --")


if(BUILD_TRELLIS)
  add_subdirectory(Trellis                                    "${CMAKE_CURRENT_BINARY_DIR}/Trellis")
endif(BUILD_TRELLIS)


if(BUILD_XOCTREE)
  add_subdirectory(xOctree                                    "${CMAKE_CURRENT_BINARY_DIR}/xOctree")
endif(BUILD_XOCTREE)



add_subdirectory(xAnalyticalSolution                        "${CMAKE_CURRENT_BINARY_DIR}/xAnalyticalSolution")
add_subdirectory(xFastMarching                      "${CMAKE_CURRENT_BINARY_DIR}/xFastMarching")
add_subdirectory(xGraph                             "${CMAKE_CURRENT_BINARY_DIR}/xGraph")
add_subdirectory(xTool                              "${CMAKE_CURRENT_BINARY_DIR}/xTool")
add_subdirectory(xTensor                            "${CMAKE_CURRENT_BINARY_DIR}/xTensor")
add_subdirectory(xGeom/xGeomTools                   "${CMAKE_CURRENT_BINARY_DIR}/xGeom/xGeomTools")
add_subdirectory(xGeom/xScanElement                   "${CMAKE_CURRENT_BINARY_DIR}/xGeom/xScanElement")
add_subdirectory(xGeom/xDistanceNearest                   "${CMAKE_CURRENT_BINARY_DIR}/xGeom/xDistanceNearest")
add_subdirectory(xMeshTool/xReader                   "${CMAKE_CURRENT_BINARY_DIR}/xMeshTool/xReader")
add_subdirectory(xMeshTool/xTemplateRefineMesh        "${CMAKE_CURRENT_BINARY_DIR}/xMeshTool/xTemplateRefineMesh")
add_subdirectory(xMeshTool/xDistMesh        "${CMAKE_CURRENT_BINARY_DIR}/xMeshTool/xDistMesh")
add_subdirectory(xMapping                            "${CMAKE_CURRENT_BINARY_DIR}/xMapping")
add_subdirectory(xQuadrature                              "${CMAKE_CURRENT_BINARY_DIR}/xQuadrature")
add_subdirectory(xLinAlg                            "${CMAKE_CURRENT_BINARY_DIR}/xLinAlg")
add_subdirectory(xInterface/AOMD                              "${CMAKE_CURRENT_BINARY_DIR}/xInterface/AOMD")
add_subdirectory(xInterface/xOctree                              "${CMAKE_CURRENT_BINARY_DIR}/xInterface/xOctree")
add_subdirectory(xInterface/xTemplateRefineMesh                              "${CMAKE_CURRENT_BINARY_DIR}/xInterface/xTemplateRefineMesh")
add_subdirectory(xInterface/Graph                              "${CMAKE_CURRENT_BINARY_DIR}/xInterface/Graph")
add_subdirectory(xFEM                              "${CMAKE_CURRENT_BINARY_DIR}/xFEM")
add_subdirectory(xExport                              "${CMAKE_CURRENT_BINARY_DIR}/xExport")
add_subdirectory(xCut/xLegacySimpleCut                              "${CMAKE_CURRENT_BINARY_DIR}/xCut/xLegacySimpleCut")
add_subdirectory(xCrack                             "${CMAKE_CURRENT_BINARY_DIR}/xCrack") 
