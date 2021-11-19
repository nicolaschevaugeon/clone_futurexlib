enable_testing()
add_test( 
  NAME xToDense 
  COMMAND ${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH}
 )
add_test( 
  NAME ndiff_xToDense 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xToDense PROPERTIES DEPENDS xToDense ) 
