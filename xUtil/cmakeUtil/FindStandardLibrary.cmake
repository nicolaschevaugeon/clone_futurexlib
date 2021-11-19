#  Find general
include(FindNoHeaderLibrary)
function(find_standard_library)

  set(libname    "${ARGV0}" )
  set(LIB_NAME   "${ARGV1}")
  set(LIB_PATH   "${ARGV2}")
  set(headerfile "${ARGV3}" )
  set(INC_NAME   "${ARGV4}")
  set(INC_PATH   "${ARGV5}")
  set(REQUIRED   "${ARGV6}")
  
  #   message("lib name = lib${libname}.*") 
  #   message("cmake variable for the lib:  ${LIB_NAME}") 
  #   message("root path to find the lib: ${${LIB_PATH}}") 
  #   message("cmake variable for the header file: ${headerfile}") 
  #   message("cmake variable for the include path: ${INC_NAME}") 
  #   message("parent path to find the header file: ${${INC_PATH}}") 
  #  message("required option: ${REQUIRED}") 

  # searching lib
  if(${LIB_NAME})

    message(STATUS "${libname} = ${${LIB_NAME}}" )

  else()

    find_no_header_library( ${libname}    ${LIB_NAME}  ${LIB_PATH}  ${REQUIRED} )

  endif()

  # searching include
  if(${INC_NAME})

    message(STATUS "includes = ${${INC_NAME}}" )

  else()

    foreach(item ${CMAKE_PREFIX_PATH} )
      get_filename_component(ROOT ${item} ABSOLUTE)
      list(APPEND POSSIBLE_PATHS "${ROOT}" )
      list(APPEND POSSIBLE_PATHS "${ROOT}/include/*")
    endforeach()
    find_path(${INC_NAME} ${headerfile}  
      HINTS ${${INC_PATH}}
      PATHS ${POSSIBLE_PATHS}
      )

    if(${INC_NAME})
      message(STATUS "Found includes: ${${INC_NAME}}" )
    else()
      message( " includes for ${headerfile} NOT FOUND" )
      if (REQUIRED)
	message( " please define ${INC_PATH} or ${INC_NAME} cmake variable with the path of include directory" )
	message(FATAL_ERROR "${headerfile} REQUIRED " )
      endif(REQUIRED)
    endif()

  endif()

  set("${ARGV1}" ${${LIB_NAME}}    CACHE FILEPATH "${libname} libraries" FORCE )
  set("${ARGV4}" ${${INC_NAME}}    CACHE PATH     "${libname} include  " FORCE )


endfunction(find_standard_library)
