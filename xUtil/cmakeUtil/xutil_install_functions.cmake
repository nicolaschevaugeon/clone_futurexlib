###########################################################################
# set LIBRARIES_INSTALL_PATH to the good value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function(set_lib_install_root)

  #define_archos_suffixe(ARCHOS)

#if("${ARGV0}" STREQUAL "")
  if( NOT ${ARGV0})
    set( LIBRARIES_INSTALL_PATH "${DEVROOT}/lib" )
  else()
    set( LIBRARIES_INSTALL_PATH "${${ARGV0}}")
  endif()

  get_filename_component(LIBRARIES_INSTALL_PATH ${LIBRARIES_INSTALL_PATH}  ABSOLUTE)
  set("${ARGV0}" ${LIBRARIES_INSTALL_PATH} PARENT_SCOPE )

endfunction(set_lib_install_root)



###########################################################################
# set LIBRARIES_INSTALL_PATH to the good value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function(install_target_in_destination)

  #define_archos_suffixe(ARCHOS)
  set_lib_install_root(LIBRARIES_INSTALL_PATH)
#  message(AUTHOR_WARNING "${ARGV0} will be installed in: ${LIBRARIES_INSTALL_PATH}/${ARCHOS}")
  install(TARGETS ${ARGV} DESTINATION ${LIBRARIES_INSTALL_PATH}/${ARCHOS})

endfunction(install_target_in_destination)
