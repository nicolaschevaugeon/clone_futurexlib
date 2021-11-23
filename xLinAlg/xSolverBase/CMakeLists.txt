cmake_minimum_required(VERSION 3.13...3.19)



message(STATUS "=============== --")
message(STATUS "entering xLinAlg/xSolverBase --")
message(STATUS "=============== --")


project(
    xSolverBase
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

find_package(BLAS         REQUIRED  )

set(XSOLVERBASE_SOURCE
# Headers   
src/xBasicIterativeSolvers.h
src/xBlasDef.h
src/xCSRMatrix.h
src/xCSRMatrix_imp.h
src/xCSRVector.h
src/xCoordinateMatrix.h
src/xDenseMatrix.h
src/xDiagonalMatrix.h
src/xDistBlockDenseMatrix.h
src/xDistBlockDenseMatrix_imp.h
src/xDistIndex.h
src/xDistVector.h
src/xDistVector_imp.h
src/xFixedCCSLowerMatrix.h
src/xFixedCRSMatrix.h
src/xGenericSparseMatrix.h
src/xGenericSparseMatrixGemvPolicy.h
src/xGenericSparseMatrixTraitPolicy.h
src/xGenericSparseMatrix_imp.h
src/xGraphMatrix.h
src/xLinearSystemSolverDiagonal.h
src/xScalapackDef.h
src/xSparseVector.h
src/xSparseVectorPackUnPack.h
src/xSparseVectorPackUnPack_imp.h
src/xSparseVector_imp.h
src/xTraitsMatrix.h
# Sources
src/xCSRMatrix.cc
src/xCSRVector.cc
src/xDenseMatrix.cc
src/xDiagonalMatrix.cc
src/xDistBlockDenseMatrix.cc
src/xDistIndex.cc
src/xFixedCRSMatrix.cc
src/xGenericSparseMatrix.cc
src/xGenericSparseMatrixTraitPolicy.cc
src/xGraphMatrix.cc
src/xLinearSystemSolverDiagonal.cc
)

add_library(xSolverBase ${BUILD_SHARED_LIBS} ${XSOLVERBASE_SOURCE})
target_include_directories(xSolverBase PUBLIC src)
add_library(xLinAlg::xSolverBase ALIAS xSolverBase)


target_link_libraries(xSolverBase PUBLIC xTool::xTool PRIVATE BLAS::BLAS)
add_dependencies(xSolverBase xTool::xTool) 



# Make alias

# Export to a Config.cmake script
export(TARGETS xSolverBase NAMESPACE xLinAlg:: FILE xSolverBaseConfig.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE xSolverBase)

# Install
#install_target_in_destination(xSolverBase)



# Tests
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
add_subdirectory(test)
#endif()
