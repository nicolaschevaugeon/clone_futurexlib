cmake_minimum_required(VERSION 3.13...3.19)
message(STATUS "=============== --")
message(STATUS "atomic test: SOLVERBASE  --")
message(STATUS "=============== --")

#define_archos_suffixe(ARCHOS)

enable_testing()

set(TEST_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/xCSRVector
    ${CMAKE_CURRENT_SOURCE_DIR}/xDenseMatrix
    ${CMAKE_CURRENT_SOURCE_DIR}/xGenericSparseMatrixCopyConstructor
    ${CMAKE_CURRENT_SOURCE_DIR}/xGenericSparseMatrixCopyTranspose
    ${CMAKE_CURRENT_SOURCE_DIR}/xGenericSparseMatrixSparseProduct
    ${CMAKE_CURRENT_SOURCE_DIR}/xGenericSparseMatrixAoperB
    ${CMAKE_CURRENT_SOURCE_DIR}/xGenericSparseMatrixGemv
    ${CMAKE_CURRENT_SOURCE_DIR}/xToDense
    ${CMAKE_CURRENT_SOURCE_DIR}/xDistIndex
    ${CMAKE_CURRENT_SOURCE_DIR}/xDistVector
    ${CMAKE_CURRENT_SOURCE_DIR}/xDistBlockDenseMatrix
    )

create_tests_from_list(LIST_TESTS ${TEST_LIST} TARGETS xLinAlg::xSolverbase xTool::xTool BLAS::BLAS)


