enable_testing()
add_test( 
  NAME xCSRVector 
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${INSTALL_PATH}/${TARGET_NAME}  data/main.dat
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_xCSRVector 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xCSRVector PROPERTIES DEPENDS xCSRVector ) 
