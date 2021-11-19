enable_testing()

#add_executable(testxOctreeGrid main.cc)
#target_link_libraries(testxOctreeGrid PRIVATE xFEM xInterfaceAOMDGeneral)


add_test(
  NAME testxOctreeGrid
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} testxOctreeGrid ${MPIEXEC_POSTFLAGS} 
  #${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

#add_test(ndiff_testxOctreeGrid ${TESTNDIFF}  )
#set_tests_properties(ndiff_testxOctreeGrid PROPERTIES DEPENDS testxOctreeGrid)
