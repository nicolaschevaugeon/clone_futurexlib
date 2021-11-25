# Find package ParMetis / Metis 
# on parallel   architecture use METIS inside PARMETIS
# on sequential architecture use stanalone METIS  
include(FindNoHeaderLibrary)
include(FindStandardLibrary)

option(USE_PARMETIS ON)

if(NOT MPI_COMPILER )
  find_package(MPI REQUIRED)
endif()
if(MPI_COMPILER)
  message(STATUS    "MPI for PARMETIS = ${MPI_COMPILER}")
  set(USE_PARMETIS ON)
else()
  message(STATUS "NO MPI for PARMETIS switch to find METIS only")
  set(USE_PARMETIS OFF)
endif() 


if(USE_PARMETIS)

  if(PARMETIS_LIBRARIES)
    message(STATUS  "parmetis + metis = ${PARMETIS_LIBRARIES}")
    message(STATUS  "                   ${METIS_LIBRARIES}")
    message(STATUS  "includes =         ${PARMETIS_INCLUDE_DIR}")
    message(STATUS  "                   ${METIS_INCLUDE_DIR}")

  else()

    find_standard_library( parmetis PARMETIS_LIBRARIES PARMETIS_LIBRARY_PATH parmetis.h PARMETIS_INCLUDE_DIR PARMETIS_INCLUDE_PATH ${METIS_FIND_REQUIRED})
    # find Metis of ParMetis
    get_filename_component(PARMETIS_LIBRARY_PATH ${PARMETIS_LIBRARIES} PATH)
    find_library(METIS_LIBRARIES metis NO_DEFAULT_PATH  HINTS ${PARMETIS_LIBRARY_PATH} ${METIS_LIBRARY_PATH}  ${METIS_FIND_REQUIRED})
    foreach(item ${CMAKE_PREFIX_PATH} )
      get_filename_component(ROOT ${item}/.. ABSOLUTE)
      list(APPEND POSSIBLE_PATHS "${ROOT}/include" )
      list(APPEND POSSIBLE_PATHS "${ROOT}/include/metis")
      list(APPEND POSSIBLE_PATHS "${ROOT}/include/*")
    endforeach()

    list(APPEND POSSIBLE_PATHS "${METIS_INCLUDE_PATH}")

    find_path(METIS_INCLUDE_DIR metis.h        PATHS ${POSSIBLE_PATHS} )

   endif()
   # check version of parmetis
   find_file(PARMETIS_HEADER parmetis.h
     HINTS ${PARMETIS_INCLUDE_DIR}
     NO_DEFAULT_PATH)
   if(PARMETIS_HEADER)
     file(READ "${PARMETIS_HEADER}" parmetisheader)
     # get version number from defines in header file
     string(REGEX REPLACE ".*#define PARMETIS_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1"
       PARMETIS_MAJOR_VERSION  "${parmetisheader}")
     string(REGEX REPLACE ".*#define PARMETIS_MINOR_VERSION[ \t]+([0-9]+).*" "\\1"
       PARMETIS_MINOR_VERSION  "${parmetisheader}")
     string(REGEX REPLACE ".*#define PARMETIS_SUBMINOR_VERSION[ \t]+([0-9]+).*" "\\1"
       PARMETIS_SUBMINOR_VERSION "${parmetisheader}")
     set(PARMETIS_VERSION ${PARMETIS_MAJOR_VERSION}.${PARMETIS_MINOR_VERSION}.${PARMETIS_SUBMINOR_VERSION})
     message(STATUS "Found Parmetis v${PARMETIS_VERSION}")
   else()
       message(STATUS "Can't find parmetis.h")
       set(PARMETIS_VERSION "0.0.0")
   endif()

else()

  if(METIS_LIBRARIES)
    message(STATUS  "metis:  ${METIS_LIBRARIES}")
  else()
    find_standard_library( metis       METIS_LIBRARIES    METIS_LIBRARY_PATH    metis.h    METIS_INCLUDE_DIR    METIS_INCLUDE_PATH  ${METIS_FIND_REQUIRED})
  endif()
  
  if(METIS_LIBRARIES)
  
endif()
endif()

  if(NOT TARGET METIS::METIS)
    add_library(METIS::METIS IMPORTED INTERFACE)
    set_property(TARGET METIS::METIS PROPERTY INTERFACE_LINK_LIBRARIES ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES})
    if(PARMETIS_INCLUDE_DIR)
     # set(ALL_INCLUDES_PARMETIS 
      set_target_properties(METIS::METIS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${PARMETIS_INCLUDE_DIR}) 
      set_target_properties(METIS::METIS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIR};${METIS_INCLUDE_DIR}")
    endif()
  endif()


