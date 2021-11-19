enable_testing()
add_test( 
  NAME test1 
  COMMAND ${TARGET_NAME}  data/main.dat
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_test1 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_test1 PROPERTIES DEPENDS test1 ) 
