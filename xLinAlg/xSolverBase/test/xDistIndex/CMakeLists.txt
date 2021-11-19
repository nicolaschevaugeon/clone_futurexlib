enable_testing()
add_test(  
  NAME xDistIndex 
  COMMAND ${MPIEXEC} -n 3 ${INSTALL_PATH}/${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xDistIndex 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xDistIndex PROPERTIES DEPENDS xDistIndex ) 
