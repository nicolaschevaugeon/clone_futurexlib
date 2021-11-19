enable_testing()
add_test(
  NAME testDataExchanger
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${INSTALL_PATH}
 )
add_test(
  NAME ndiff_testDataExchanger
  COMMAND ${TESTNDIFF}
  WORKING_DIRECTORY ${INSTALL_PATH}
  )
set_tests_properties(ndiff_testDataExchanger PROPERTIES DEPENDS testDataExchanger ) 

