enable_testing()
add_test(
  NAME  mpiInitThreadMultiple
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${INSTALL_PATH}
  )
# For MVAPICH2-2.0  need to set MV2_ENABLE_AFFINITY to 0 to get MPI_THREAD_MULTIPLE available
# For other MPI library it normaly won't have any influence
set_tests_properties(mpiInitThreadMultiple PROPERTIES ENVIRONMENT "MV2_ENABLE_AFFINITY=0")
add_test(
  NAME ndiff_mpiInitThreadMultiple
  COMMAND ${TESTNDIFF}
  WORKING_DIRECTORY ${INSTALL_PATH}
  )
set_tests_properties(ndiff_mpiInitThreadMultiple PROPERTIES DEPENDS mpiInitThreadMultiple ) 
