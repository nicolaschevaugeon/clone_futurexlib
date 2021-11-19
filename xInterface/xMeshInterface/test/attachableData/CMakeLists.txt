enable_testing()
add_test(  
  NAME atomic_attachableData 
  COMMAND  ${TARGET_NAME}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
)
add_test( 
  NAME ndiff_atomic_attachableData
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_atomic_attachableData PROPERTIES DEPENDS atomic_attachableData ) 
