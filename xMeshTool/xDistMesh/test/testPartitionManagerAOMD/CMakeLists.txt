add_test( 
  NAME testPartitionManagerAOMD 
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}  ${INSTALL_PATH}/data/test.msh
  WORKING_DIRECTORY ${INSTALL_PATH}
)
add_test( 
  NAME ndiff_testPartitionManagerAOMD 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_testPartitionManagerAOMD PROPERTIES DEPENDS testPartitionManagerAOMD ) 

