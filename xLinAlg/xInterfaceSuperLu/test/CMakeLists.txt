message(STATUS       "----------------------------------")
message(STATUS "atomic test: INTERFACE SUPERLU  --")
message(STATUS       "----------------------------------")

enable_testing()

create_tests()

set(TEST_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/simple
    )

create_tests_from_list(LIST_TESTS ${TEST_LIST} TARGETS xLinAlg::xSolverbase )
