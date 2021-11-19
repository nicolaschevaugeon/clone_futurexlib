cmake_minimum_required(VERSION 3.13..3.19)

message(STATUS       "---------------------------- --")
message(STATUS "atomic test: INTERFACE LAPACK  --")
message(STATUS       "---------------------------- --")


enable_testing()


create_tests(TARGETS xLinAlg::xInterfaceLapack)




