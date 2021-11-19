enable_testing()
add_test(
  NAME xDeltaTime
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
add_test(
  NAME ndiff_xDeltaTime
  COMMAND  ${TESTNDIFF}
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xDeltaTime PROPERTIES DEPENDS xDeltaTime ) 

