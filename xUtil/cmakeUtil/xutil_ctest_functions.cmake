set(XFILES_CMAKE_UTIL_PATH "${DEVROOT}/Xfiles/xUtil/cmakeUtil"
    CACHE PATH
    "path to Xfiles/Util/cmakeUtil")

###########################################################################
# create all executable-targets for test/* 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function(create_tests_from_list)

  # Use cmake_parse_arguments...
  set(prefix ARG)
  #set(singleValues TARGETS)
  set(multiValues LIST_TESTS TARGETS)
  
  cmake_parse_arguments(${prefix}
                        ""
                        ""
                        "${multiValues}"
                        ${ARGN})
                        

  if("${EXECUTABLES_INSTALL_PATH}" STREQUAL "")
    set(EXECUTABLES_INSTALL_PATH ${CMAKE_CURRENT_BINARY_DIR})
  else()
    set(OUT_OF_BUILD_INSTALL "TRUE")
  endif()

  message(STATUS "installation path: ${EXECUTABLES_INSTALL_PATH}")
  set(TESTNDIFF "${XFILES_CMAKE_UTIL_PATH}/test_ndiff.sh")

  message(STATUS "ARG_LIST_TESTS ${ARG_LIST_TESTS}")
    
  foreach(item ${ARG_LIST_TESTS})
    set(test ${item})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test} OR EXISTS ${test})
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test} )
    set(test ${CMAKE_CURRENT_SOURCE_DIR}/${item})
      endif()
    else()
      get_filename_component(tmp ${test} NAME )
      message(FATAL_ERROR "Directory ${CMAKE_CURRENT_SOURCE_DIR}/${tmp} does not exist!")
    endif()

    if(EXISTS ${test}/CMakeLists.txt )
      get_filename_component(testname ${test} NAME) 

      if (OUT_OF_BUILD_INSTALL)
        set(TARGET_NAME ${testname}_${ARCHOS})
        set(INSTALL_PATH ${EXECUTABLES_INSTALL_PATH}/${testname})
      else()
        set(TARGET_NAME ${testname})
        set(INSTALL_PATH ${CMAKE_CURRENT_BINARY_DIR}/${testname})
      endif()
      if(EXISTS ${test}/main.cc OR EXISTS ${test}/main.cpp)
        message("   -- create test: ${testname} in ${INSTALL_PATH} ")
        if(EXISTS ${test}/main.cc)
          add_executable("${TARGET_NAME}" ${CMAKE_CURRENT_SOURCE_DIR}/${testname}/main.cc )
        endif()
        if(EXISTS ${test}/main.cpp )
          add_executable("${TARGET_NAME}" ${CMAKE_CURRENT_SOURCE_DIR}/${testname}/main.cpp )
        endif()
        target_link_libraries("${TARGET_NAME}" ${EXTERNAL_LIBRARIES} ${ARG_TARGETS})
        message(STATUS "ARG_TARGET is ${ARG_TARGETS}")
        set_target_properties("${TARGET_NAME}"
          PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY ${INSTALL_PATH}
          COMPILE_FLAGS "-Wno-deprecated"
          )
        add_custom_command(TARGET "${TARGET_NAME}"
              COMMAND  bash "${XFILES_CMAKE_UTIL_PATH}/do_link.sh" ${CMAKE_CURRENT_SOURCE_DIR}/$        {testname}/data       ${EXECUTABLES_INSTALL_PATH}/${testname}/data)

        add_custom_command(TARGET "${TARGET_NAME}"
          COMMAND  bash "${XFILES_CMAKE_UTIL_PATH}/do_link.sh" ${CMAKE_CURRENT_SOURCE_DIR}/${testname}/reference  ${EXECUTABLES_INSTALL_PATH}/${testname}/reference)
    
      endif()
      add_subdirectory(${testname})
      enable_testing()      
    endif()
  endforeach()
endfunction(create_tests_from_list)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





###############################################################################################
# create all executable-targets for test/xlinalg (same src but not same solver) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
function(create_solver_tests)

  # Use cmake_parse_arguments...
  set(prefix ARG)
  #set(singleValues TARGETS)
  set(multiValues LIST_TESTS TARGETS)
  
  cmake_parse_arguments(${prefix}
                        ""
                        ""
                        "${multiValues}"
                        ${ARGN})

  message(STATUS "ARG_LIST_TESTS ${ARG_LIST_TESTS}")
  

  get_filename_component(solvername ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  if("${EXECUTABLES_INSTALL_PATH}" STREQUAL "")
    set(EXECUTABLES_INSTALL_PATH ${CMAKE_CURRENT_BINARY_DIR})
  else()
    set(OUT_OF_BUILD_INSTALL "TRUE")
  endif()

  message(STATUS "installation path: ${EXECUTABLES_INSTALL_PATH}")
  #message_verbose(text)
  #define_archos_suffixe(ARCHOS)
  set(TESTNDIFF "${XFILES_CMAKE_UTIL_PATH}/test_ndiff.sh")

#  message("ARGV: ${ARGV}")
  #foreach(test ${ARGV} )
  foreach(test ${ARG_LIST_TESTS})   
    get_filename_component(testname ${test} NAME) 
    set(testname ${testname}_${solvername})
    if (OUT_OF_BUILD_INSTALL)
      set(TARGET_NAME ${testname}_${ARCHOS})
      set(INSTALL_PATH ${EXECUTABLES_INSTALL_PATH}/${testname})
    else()
      set(TARGET_NAME ${testname})
      set(INSTALL_PATH ${CMAKE_CURRENT_BINARY_DIR}/${testname})
    endif()

    message("   -- create test : ${testname} ")
    add_executable("${TARGET_NAME}"  ${test}/main.cc )
    #include_directories( ${EXTERNAL_INCLUDES} )
    #target_link_libraries("${TARGET_NAME}" ${EXTERNAL_LIBRARIES})
    target_link_libraries("${TARGET_NAME}" ${EXTERNAL_LIBRARIES} ${ARG_TARGETS})

    set_target_properties("${TARGET_NAME}"
      PROPERTIES  
      RUNTIME_OUTPUT_DIRECTORY ${INSTALL_PATH}
      COMPILE_FLAGS "-Wno-deprecated"
      )

    add_custom_command(TARGET "${TARGET_NAME}"
      COMMAND  "${XFILES_CMAKE_UTIL_PATH}/do_link.sh" ${test}/data  ${INSTALL_PATH}/data)

    if(EXISTS ${test}/reference )
      add_custom_command(TARGET "${TARGET_NAME}"
        COMMAND  "${XFILES_CMAKE_UTIL_PATH}/do_link.sh" ${test}/reference  ${INSTALL_PATH}/reference)
    else()
      add_custom_command(TARGET "${TARGET_NAME}"
        COMMAND  "${XFILES_CMAKE_UTIL_PATH}/do_link.sh" ${CMAKE_CURRENT_SOURCE_DIR}/${testname}/reference  ${INSTALL_PATH}/reference)
    endif()

    add_subdirectory(${testname})

  endforeach()


  
  
endfunction(create_solver_tests)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###########################################################################
# create all executable-targets for test/* 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function(create_tests)

  # Use cmake_parse_arguments...
  set(prefix ARG)
  #set(singleValues TARGETS)
  set(multiValues TARGETS)
  
  cmake_parse_arguments(${prefix}
                        ""
                        ""
                        "${multiValues}"
                        ${ARGN})

  message(STATUS "ARG_TARGETS() ${ARG_TARGETS}")

  execute_process( COMMAND chmod u+x ${XFILES_CMAKE_UTIL_PATH}/do_link.sh )
  execute_process( COMMAND chmod u+x ${XFILES_CMAKE_UTIL_PATH}/test_ndiff.sh )
  execute_process( COMMAND chmod u+x ${XFILES_CMAKE_UTIL_PATH}/purge.sh )

  file(GLOB  TESTS * )
  set(LIST ${TESTS})
  remove(LIST 
    ${CMAKE_CURRENT_SOURCE_DIR}/NOTICE 
    ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE 
    ${CMAKE_CURRENT_SOURCE_DIR}/Testing 
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake* 
    ${CMAKE_CURRENT_SOURCE_DIR}/Make*
    ${CMAKE_CURRENT_SOURCE_DIR}/.svn
    ${CMAKE_CURRENT_SOURCE_DIR}/.hg
    ${CMAKE_CURRENT_SOURCE_DIR}/.hgignore
    )
  get_filename_component(MODULE ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)
  get_filename_component(SUBMODULE ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
  get_filename_component(MODULE ${MODULE} NAME)
  get_filename_component(SUBMODULE ${SUBMODULE} NAME)
  message(STATUS "${MODULE}/${SUBMODULE}:")
  #message_verbose(text)

  message (STATUS "ATT ${ARG_TARGETS}")
  message (STATUS "ATL ${LIST}")
  create_tests_from_list(LIST_TESTS ${LIST} TARGETS ${ARG_TARGETS})

endfunction(create_tests)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
