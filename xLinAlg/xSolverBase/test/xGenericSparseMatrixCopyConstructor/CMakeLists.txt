enable_testing()
add_test(  
  NAME xGenericSparseMatrixCopyConstructor 
  COMMAND ${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xGenericSparseMatrixCopyConstructor 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xGenericSparseMatrixCopyConstructor PROPERTIES DEPENDS xGenericSparseMatrixCopyConstructor ) 
