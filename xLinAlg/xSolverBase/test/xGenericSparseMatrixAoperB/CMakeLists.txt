enable_testing()
add_test( 
  NAME xGenericSparseMatrixAoperB 
  COMMAND ${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xGenericSparseMatrixAoperB 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xGenericSparseMatrixAoperB PROPERTIES DEPENDS xGenericSparseMatrixAoperB ) 
