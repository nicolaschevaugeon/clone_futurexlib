enable_testing()
add_test(  
  NAME xGenericSparseMatrixCopyTranspose 
  COMMAND ${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
 )
add_test( 
  NAME ndiff_xGenericSparseMatrixCopyTranspose 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xGenericSparseMatrixCopyTranspose PROPERTIES DEPENDS xGenericSparseMatrixCopyTranspose ) 
