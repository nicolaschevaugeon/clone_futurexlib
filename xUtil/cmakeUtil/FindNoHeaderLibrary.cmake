function(find_no_header_library)

  set(libname    "${ARGV0}" )
  set(LIB_NAME   "${ARGV1}")
  set(LIB_PATH   "${ARGV2}")
  set(REQUIRED   "${ARGV3}")
  
  #  message("lib name = lib${libname}.*") 
  #  message("cmake variable for the lib: ${LIB_NAME}") 
  #  message("cmake variable for the lib path: ${LIB_PATH}") 
  # message("required option: ${REQUIRED}") 

  if(${LIB_NAME})

    set(text "${libname}: ${${LIB_NAME}}" )
    message_verbose(text)

  else(${LIB_NAME})

    #set(text "looking for ${libname} ..." )
    # message_verbose(text)
    foreach(item ${CMAKE_PREFIX_PATH} )
      get_filename_component(ROOT ${item} ABSOLUTE)
      list(APPEND POSSIBLE_PATHS "${ROOT}" )
      list(APPEND POSSIBLE_PATHS "${ROOT}/*" )
    endforeach()


    
    find_library(${LIB_NAME}  ${libname} 
      NAME  lib${libname}.a 
            lib${libname}.so 
      HINTS ${${LIB_PATH}}
      PATHS  
      ${POSSIBLE_PATHS}
      )
    
    if(${LIB_NAME})
      message(STATUS "Found ${libname}: ${${LIB_NAME}}" )
    else()
      message( " ${libname} NOT FOUND" )
      message( " please define ${LIB_PATH} cmake variable with the path of the directory" )
      message( " or define ${LIB_NAME} with the path of the library" )

      if("${REQUIRED}" STREQUAL "1")
        message(FATAL_ERROR "${libname} REQUIRED" )
      endif()
      
    endif()

  endif(${LIB_NAME})

  #  set("${ARGV1}" ${${LIB_NAME}}   PARENT_SCOPE)
  set("${ARGV1}" ${${LIB_NAME}}    CACHE FILEPATH "${libname} libraries" FORCE )

endfunction(find_no_header_library)
