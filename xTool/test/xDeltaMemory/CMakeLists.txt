enable_testing()
add_test(
  NAME xDeltaMemory
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
add_test(
  NAME ndiff_xDeltaMemory
  COMMAND  ${TESTNDIFF}
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_xDeltaMemory PROPERTIES DEPENDS xDeltaMemory ) 

