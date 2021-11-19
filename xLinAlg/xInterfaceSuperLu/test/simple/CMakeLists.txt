enable_testing()
add_test(  
  NAME simpleSuperLU 
  COMMAND ${TARGET_NAME}  data/main.dat
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test(ndiff_simpleSuperLU ${TESTNDIFF}  )
set_tests_properties(ndiff_simpleSuperLU PROPERTIES DEPENDS simpleSuperLU ) 
