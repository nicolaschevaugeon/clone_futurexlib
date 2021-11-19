enable_testing()
add_test(   
  NAME xDistBlockDenseMatrix 
  COMMAND ${MPIEXEC} -n 4 ${INSTALL_PATH}/${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xDistBlockDenseMatrix 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xDistBlockDenseMatrix PROPERTIES DEPENDS xDistBlockDenseMatrix ) 
