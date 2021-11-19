#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2021 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# --
#[=======================================================================[.rst:
FindMUMPS
-----------

Find mumps libraries and headers

Usage :
 
find_package(MUMPS REQUIRED  COMPONENTS <name>)
target_link_libraries(yourlib PRIVATE MUMPS::MUMPS)

with name among :
- SEQ for sequential version of mumps
- PAR for paralle version

If COMPONENTS is not set, SEQ version is searched.
Warning : do not use several components in your build.

Notice that some extra libs are also searched for : ptscotch scotch metis pord.

Set MUMPS_ROOT=<where mumps is installed>
if it's not in a "classic" place or if you want a specific version

header : dmumps_c.h
libs : see MUMPS_NAME var below.

#]=======================================================================]

include(FindPackageHandleStandardArgs)



 #Try to provide some hints for a find_package call.
#
# Usage:
#
# set_find_package_hints(NAME <name> MODULE <mod>)
#
#
# Result : set (parent scope) _<NAME>_SEARCH_OPTS and _<NAME>_INC_SEARCH_OPTS
# that can be used in find_path (INC_SEARCH) and find_library calls.
#
# These variables are filled either with <NAME>_ROOT value if set
# or using pkg-config information, if available.
#
# See examples of use in FindCPPUNIT.cmake or FindSuperLU.cmake.
#
function(set_find_package_hints)
  set(oneValueArgs NAME MODULE)

  cmake_parse_arguments(pkg "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(${pkg_NAME}_ROOT)
    set(_${pkg_NAME}_SEARCH_OPTS
      HINTS ${${pkg_NAME}_ROOT} NO_DEFAULT_PATH PARENT_SCOPE)
    set(_${pkg_NAME}_INC_SEARCH_OPTS
      HINTS ${${pkg_NAME}_ROOT} NO_DEFAULT_PATH PARENT_SCOPE)
  else()
    # Try pkgconfig
    find_package(PkgConfig QUIET)
    pkg_check_modules(PKGC_${pkg_NAME} ${pkg_MODULE} QUIET)
    if(PKGC_${pkg_NAME}_FOUND)
      set(_${pkg_NAME}_INC_SEARCH_OPTS "HINTS ${PKGC_${pkg_NAME}_INCLUDE_DIRS}"
        PARENT_SCOPE)
    endif()
    set(_${pkg_NAME}_SEARCH_OPTS
      HINTS ${PKGC_${pkg_NAME}_LIBRARY_DIRS} ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
      PARENT_SCOPE)
  endif()

endfunction()


list(LENGTH MUMPS_FIND_COMPONENTS nb_comp)
if(nb_comp GREATER 1)
  message(FATAL_ERROR "Try to find several MUMPS conflicting components : ${MUMPS_FIND_COMPONENTS}.")
endif()
if(NOT MUMPS_FIND_COMPONENTS)
  set(MUMPS_FIND_COMPONENTS SEQ)
endif()

if(${MUMPS_FIND_COMPONENTS} STREQUAL SEQ)
  set(__MUMPS_NAMES dmumps_seq dmumps esmumps)
  set(__SUFFIX seq)
  # on debian systems one may have mumps+[pt]scotch packages also
  set(__MUMPS_COMMON_NAMES mumps_common_seq mumps_common esmumps)
elseif(${MUMPS_FIND_COMPONENTS} STREQUAL PAR)
  set(__MUMPS_NAMES dmumps_ptscotch dmumps_scotch dmumps esmumps)
  set(__MUMPS_COMMON_NAMES mumps_common_ptscotch mumps_common_scotch mumps_common esmumps)
else()
  message(FATAL_ERROR "Unknown MUMPS component ${MUMPS_FIND_COMPONENTS}. Search failed.")
endif()



if(NOT MUMPS_ROOT)
  set(MUMPS_ROOT $ENV{MUMPS_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME MUMPS MODULE mumps)

find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h
  PATH_SUFFIXES include mumps mumps${__SUFFIX}
  ${_MUMPS_INC_SEARCH_OPTS}
  )

if(NOT MUMPS_LIBRARIES)
  find_library(MUMPS_LIBRARY NAMES ${__MUMPS_NAMES}
    ${_MUMPS_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
  if(MUMPS_LIBRARY)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_LIBRARY})
  endif()
  find_library(MUMPS_COMMON_LIBRARY NAMES ${__MUMPS_COMMON_NAMES}
    ${_MUMPS_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
  if(MUMPS_COMMON_LIBRARY)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_COMMON_LIBRARY})
  endif()
endif()

# -- Library setup --
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)

# Extras
set(extras_libs ptscotch scotch scotcherr metis pord${__SUFFIX} scalapack esmumps)
foreach(extra IN LISTS extras_libs)
  find_library(MUMPS_${extra}_LIBRARY NAMES ${extra}
    ${_MUMPS_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
  if(MUMPS_${extra}_LIBRARY)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_${extra}_LIBRARY})
    message(STATUS "Found ${MUMPS_${extra}_LIBRARY} --> ${MUMPS_LIBRARIES}")
  endif()
    
endforeach()



if(MUMPS_FOUND)
  
  if(NOT TARGET MUMPS::MUMPS)
    add_library(MUMPS::MUMPS IMPORTED INTERFACE)
    set_property(TARGET MUMPS::MUMPS PROPERTY INTERFACE_LINK_LIBRARIES ${MUMPS_LIBRARIES})
    if(MUMPS_INCLUDE_DIR)
      set_target_properties(MUMPS::MUMPS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}")
    endif()
  endif()
endif()

