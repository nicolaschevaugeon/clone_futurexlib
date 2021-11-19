enable_testing()
add_test( 
  NAME xGenericSparseMatrixSparseProduct 
  COMMAND ${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xGenericSparseMatrixSparseProduct 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xGenericSparseMatrixSparseProduct PROPERTIES DEPENDS xGenericSparseMatrixSparseProduct ) 
