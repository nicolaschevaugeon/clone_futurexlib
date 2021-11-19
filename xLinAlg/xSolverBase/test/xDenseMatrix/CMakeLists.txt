enable_testing()
add_test( 
  NAME xDenseMatrix 
  COMMAND ${TARGET_NAME} data/main.dat
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xDenseMatrix 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xDenseMatrix PROPERTIES DEPENDS xDenseMatrix ) 
