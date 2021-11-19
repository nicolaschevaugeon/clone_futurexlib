enable_testing()
add_test(
  NAME testxRegularGrid
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

add_test(ndiff_testxRegularGrid ${TESTNDIFF}  )
set_tests_properties(ndiff_testxRegularGrid PROPERTIES DEPENDS testxRegularGrid)
