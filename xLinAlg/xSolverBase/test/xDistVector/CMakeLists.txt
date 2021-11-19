enable_testing()
add_test(   
  NAME xDistVector 
  COMMAND ${MPIEXEC} -n 3 ${INSTALL_PATH}/${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xDistVector 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xDistVector PROPERTIES DEPENDS xDistVector ) 
