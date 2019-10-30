#.rst:
# FindLAPACKE
# -------------
#
# Find LAPACKE include dirs and libraries
#
# Use this module by invoking find_package with the form::
#
#   find_package(LAPACKE
#     [REQUIRED]             # Fail with error if LAPACKE is not found
#     [COMPONENTS <libs>...] # List of libraries to look for
#     )
#
# Valid names for COMPONENTS libraries are::
#
#   ALL                      - Find all libraries
#   LAPACKE_H                - Find the lapacke.h header file
#   LAPACKE                  - Find a LAPACKE library
#   LAPACK                   - Find a LAPACK library
#   CBLAS                    - Find a CBLAS library
#   BLAS                     - Find a BLAS library
#
#  Not specifying COMPONENTS is identical to choosing ALL
#
# This module defines::
#
#   LAPACKE_FOUND            - True if headers and requested libraries were found
#   LAPACKE_INCLUDE_DIRS     - LAPACKE include directories
#   LAPACKE_LIBRARIES        - LAPACKE component libraries to be linked
#
#
# This module reads hints about search locations from variables
# (either CMake variables or environment variables)::
#
#   LAPACKE_ROOT             - Preferred installation prefix for LAPACKE
#   LAPACKE_DIR              - Preferred installation prefix for LAPACKE
#
#
# The following :prop_tgt:`IMPORTED` targets are also defined::
#
#   LAPACKE::LAPACKE         - Imported target for the LAPACKE library
#   LAPACKE::LAPACK          - Imported target for the LAPACK library
#   LAPACKE::CBLAS           - Imported target for the CBLAS library
#   LAPACKE::BLAS            - Imported target for the BLAS library
#

# ==============================================================================
# Copyright (c) 2019      The University of Tennessee and The University
#                         of Tennessee Research Foundation.  All rights
#                         reserved.
#
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
# Author: Aurelien Bouteiller <bouteill@icl.utk.edu>
# ==============================================================================

# ==============================================================================
# Copyright 2018 Damien Nguyen <damien.nguyen@alumni.epfl.ch>
#
# Distributed under the OSI-approved BSD License (the "License")
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# ==============================================================================

set(LAPACKE_SEARCH_PATHS
  ${LAPACKE_ROOT}
  $ENV{LAPACKE_ROOT}
  ${BLAS_LIBRARY}
  ${LAPACKE_DIR}
  $ENV{LAPACKE_DIR}
  ${CMAKE_PREFIX_PATH}
  $ENV{CMAKE_PREFIX_PATH}
  /usr
  /usr/local/
  /usr/local/opt # homebrew on mac
  /opt
  /opt/local
  /opt/LAPACKE
  )

set(LIB_PATH_SUFFIXES
  lib64
  lib
  lib/x86_64-linux-gnu
  lib32
  )

set(INC_PATH_SUFFIXES
  include
  include/lapack
  include/lapacke/
  lapack/include
  lapacke/include
  CBLAS/include
  LAPACKE/include
  )

if(APPLE)
  list(APPEND LIB_PATH_SUFFIXES lapack/lib openblas/lib)
elseif(WIN32)
  list(APPEND LAPACKE_SEARCH_PATHS "C:/Program Files (x86)/LAPACK")
  list(APPEND LAPACKE_SEARCH_PATHS "C:/Program Files/LAPACK")
endif()

# ==============================================================================
# Prepare some helper variables

set(LAPACKE_INCLUDE_DIRS)
set(LAPACKE_LIBRARIES)
set(LAPACKE_REQUIRED_VARS)
set(LAPACKE_FIND_ALL_COMPONENTS 0)

# ==============================================================================

macro(_find_library_with_header component incname)
  find_library(LAPACKE_${component}_LIB
    NAMES ${ARGN}
    NAMES_PER_DIR
    PATHS ${LAPACKE_SEARCH_PATHS}
    PATH_SUFFIXES ${LIB_PATH_SUFFIXES})
  if(LAPACKE_${component}_LIB)    
    set(LAPACKE_${component}_LIB_FOUND 1)
  endif()
  list(APPEND LAPACKE_REQUIRED_VARS "LAPACKE_${component}_LIB")

  # If necessary, look for the header file as well
  if(NOT "${incname}" STREQUAL "")
    find_path(LAPACKE_${component}_INCLUDE_DIR
      NAMES ${incname}
      PATHS ${LAPACKE_SEARCH_PATHS}
      PATH_SUFFIXES ${INC_PATH_SUFFIXES})
    list(APPEND LAPACKE_REQUIRED_VARS "LAPACKE_${component}_INCLUDE_DIR")
    if(LAPACKE_${component}_LIB)
      set(LAPACKE_${component}_INC_FOUND 1)
    endif()
  else()
    set(LAPACKE_${component}_INC_FOUND 1)
  endif()

  if(LAPACKE_${component}_LIB_FOUND AND LAPACKE_${component}_INC_FOUND)
    set(LAPACKE_${component}_FOUND 1)
  else()
    set(LAPACKE_${component}_FOUND 0)
  endif()
endmacro()

# ------------------------------------------------------------------------------

if(NOT LAPACKE_FIND_COMPONENTS OR LAPACKE_FIND_COMPONENTS STREQUAL "ALL")
  set(LAPACKE_FIND_ALL_COMPONENTS 1)
  set(LAPACKE_FIND_COMPONENTS "LAPACKE;LAPACK;CBLAS;BLAS")
endif(NOT LAPACKE_FIND_COMPONENTS OR LAPACKE_FIND_COMPONENTS STREQUAL "ALL")

# Make sure that all components are in capitals
set(_tmp_component_list)
foreach(_comp ${LAPACKE_FIND_COMPONENTS})
  string(TOUPPER ${_comp} _comp)
  list(APPEND _tmp_component_list ${_comp})
endforeach()
set(LAPACKE_FIND_COMPONENTS ${_tmp_component_list})
set(_tmp_component_list)
list(SORT LAPACKE_FIND_COMPONENTS) # Find BLAS, then CBLAS, then LAPACK, then LAPACKE

if("${CMAKE_C_COMPILER_ID}" MATCHES ".*Clang.*" OR
   "${CMAKE_C_COMPILER_ID}" MATCHES ".*GNU.*" OR
   "${CMAKE_C_COMPILER_ID}" MATCHES ".*Intel.*"
    ) #NOT MSVC
  set(MATH_LIB "m")
 endif()

# First try the FindBLAS package
find_package(BLAS)
if(BLAS_FOUND)
  if("${BLA_VENDOR}" STREQUAL "IBMESSL")
    # Look for <essl.h>
    string(REGEX REPLACE "/lib6?4?/?[^/;]*" "" BLAS_essl_DIRS ${BLAS_LIBRARIES})
    find_path(BLAS_essl_INCLUDE_DIR
      NAMES essl.h
      PATHS ${BLAS_INCLUDE_DIRS} ${BLAS_essl_LIBRARY} ${BLAS_LIBRARY} ${BLAS_essl_DIRS}
      PATH_SUFFIXES ${INC_PATH_SUFFIXES})
    if(BLAS_essl_INCLUDE_DIR)
      list(APPEND BLAS_INCLUDE_DIRS ${BLAS_essl_INCLUDE_DIR})
      list(REMOVE_DUPLICATES BLAS_INCLUDE_DIRS)
    endif()
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${MATH_LIB})
    set(CMAKE_REQUIRED_INCLUDES ${BLAS_INCLUDE_DIRS})
    check_symbol_exists(zgemm "essl.h" FOUND_ESSL_H)
    cmake_pop_check_state()
    if(NOT FOUND_ESSL_H)
      message(WARNING "BLA_VENDOR=IBMESSL but 'essl.h' could not be found. Set BLAS_INCLUDE_DIRS by hand, or select another vendor in BLA_VENDOR")
    endif()
  endif()
  cmake_push_check_state()
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${MATH_LIB})
  set(CMAKE_REQUIRED_INCLUDES ${BLAS_INCLUDE_DIRS})
  if(CMAKE_Fortran_COMPILER_WORKS)
    check_fortran_function_exists(zgeqrf BLAS_HAS_LAPACK)
  else()
    check_c_source_compiles("int main(void) { zgeqrf_(); return 0; }" BLAS_HAS_LAPACK)
  endif()
  check_c_source_compiles("int main(void) { cblas_zgemm(); return 0; }" BLAS_HAS_CBLAS)
  check_c_source_compiles("int main(void) { LAPACKE_zgeqrf(); return 0; }" BLAS_HAS_LAPACKE)
  cmake_pop_check_state()
endif(BLAS_FOUND)

# Ok, now look for a ref-lapack
foreach(_comp ${LAPACKE_FIND_COMPONENTS})
  if(_comp STREQUAL "LAPACKE")
    if("${BLA_VENDOR}" STREQUAL "IBMESSL" OR NOT BLAS_HAS_LAPACKE)
      # LAPACKE is not in BLAS, or we have IBMESSL
      # As IBMESSL is incomplete, we complete it with ref lapack
      _find_library_with_header(${_comp} lapacke.h lapacke liblapacke)
    else()
      set(LAPACKE_LAPACKE_FOUND 1)
      set(LAPACKE_LAPACKE_LIB_FOUND 1)
    endif()
  elseif(_comp STREQUAL "LAPACKE_H")
    find_path(LAPACKE_${_comp}_INCLUDE_DIR
      NAMES lapacke.h
      PATHS ${LAPACKE_SEARCH_PATHS}
      PATH_SUFFIXES include lapack/include)
    list(APPEND LAPACKE_REQUIRED_VARS "LAPACKE_${_comp}_INCLUDE_DIR")
    if(LAPACKE_${_comp}_LIB)
      set(LAPACKE_${_comp}_INC_FOUND 1)
    endif()
  elseif(_comp STREQUAL "LAPACK")
    if("${BLA_VENDOR}" STREQUAL "IBMESSL" OR NOT BLAS_HAS_LAPACK)
      # LAPACK is not in BLAS, or we have IBMESSL
      # As IBMESSL is incomplete, we complete it with ref lapack
      _find_library_with_header(${_comp} "" lapack liblapack)
    else()
      set(LAPACKE_LAPACK_FOUND 1)
      set(LAPACKE_LAPACK_LIB_FOUND 1)
    endif()
  elseif(_comp STREQUAL "CBLAS")
    if(NOT BLAS_HAS_CBLAS)
      # CBLAS is not in BLAS, or we have IBMESSL
      # As IBMESSL is incomplete, we complete it with ref lapack
      _find_library_with_header(${_comp} cblas.h cblas libcblas)
    else()
      set(LAPACKE_CBLAS_FOUND 1)
      set(LAPACKE_CBLAS_INC_FOUND 1)
      set(LAPACKE_CBLAS_LIB_FOUND 1)
    endif()
  elseif(_comp STREQUAL "BLAS")
    if(NOT BLAS_FOUND)
      _find_library_with_header(${_comp} "" blas refblas)
      set(BLA_VENDOR CACHE "Generic")
    else()
      set(LAPACKE_BLAS_FOUND 1)
      set(LAPACKE_BLAS_LIB_FOUND 1)
      if(NOT "${BLAS_LIBRARIES}" STREQUAL "")
        set(LAPACKE_BLAS_LIB "${BLAS_LIBRARIES}")
        list(APPEND LAPACKE_REQUIRED_VARS "LAPACKE_BLAS_LIB")
      else()
        list(APPEND LAPACKE_REQUIRED_VARS "BLAS_FOUND")
      endif()
      if(NOT "${BLAS_INCLUDE_DIRS}" STREQUAL "")
        set(LAPACKE_INCLUDE_DIR "${BLAS_INCLUDE_DIRS}")
        list(APPEND LAPACKE_REQUIRED_VARS "LAPACKE_INCLUDE_DIR")
      endif()
    endif()
  else()
    message(FATAL_ERROR "Unknown component: ${_comp}")
  endif()
  mark_as_advanced(
    LAPACKE_${_comp}_LIB
    LAPACKE_${_comp}_INCLUDE_DIR)
endforeach()

# ==============================================================================

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE
  FOUND_VAR LAPACKE_FOUND
  REQUIRED_VARS ${LAPACKE_REQUIRED_VARS}
  HANDLE_COMPONENTS)

# ==============================================================================

if(LAPACKE_FOUND)
  if("${BLA_VENDOR}" STREQUAL "IBMESSL")
    # Force using ESSL first, fallback to ref-lapack if function is not found
    list(APPEND LAPACKE_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
    list(APPEND LAPACKE_LIBRARIES ${BLAS_LIBRARIES})
  endif()

  list(REVERSE LAPACKE_FIND_COMPONENTS) # For static link, have blas.a come last
  foreach(_comp ${LAPACKE_FIND_COMPONENTS})
    list(APPEND LAPACKE_INCLUDE_DIRS ${LAPACKE_${_comp}_INCLUDE_DIR})
    list(APPEND LAPACKE_LIBRARIES ${LAPACKE_${_comp}_LIB})
  endforeach()
  list(APPEND LAPACKE_LIBRARIES ${MATH_LIB})

#  if("${BLA_VENDOR}" STREQUAL "IBMESSL")
#    # ref-lapack should be compiled to use ESSL BLAS as well; add a duplicate at
#    # the end of the link chain
#    list(APPEND LAPACKE_LIBRARIES ${BLAS_LIBRARIES})
#  endif()

  if(NOT "${LAPACKE_INCLUDE_DIRS}" STREQUAL "")
    list(REMOVE_DUPLICATES LAPACKE_INCLUDE_DIRS)
  endif()

  # ----------------------------------------------------------------------------

  # Inspired by FindBoost.cmake
  foreach(_comp ${LAPACKE_FIND_COMPONENTS})
    if(LAPACKE_${_comp}_FOUND AND NOT TARGET LAPACKE::${_comp})
        #if("${BLA_VENDOR}" STREQUAL "IBMESSL")
      if(BLAS_HAS_${_comp} OR "${_comp}" STREQUAL "BLAS")
        list(LENGTH BLAS_LIBRARIES _len)
        if(_len GREATER 1)
          list(SUBLIST BLAS_LIBRARIES 1 -1 _fallback_${_comp})
        else()
          set(_fallback_${_comp})
        endif()
        if(_len GREATER 0)
          list(APPEND _fallback_${_comp} "${LAPACKE_${_comp}_LIB}")
          list(GET BLAS_LIBRARIES 0 LAPACKE_${_comp}_LIB)
        endif()
      endif()
      get_filename_component(LIB_EXT "${LAPACKE_${_comp}_LIB}" EXT)
      if(LIB_EXT STREQUAL "" OR LIB_EXT STREQUAL ".framework")
        set(LIB_TYPE INTERFACE)
      elseif(LIB_EXT STREQUAL ".a" OR LIB_EXT STREQUAL ".lib")
        set(LIB_TYPE STATIC)
      else()
        set(LIB_TYPE SHARED)
      endif()
      add_library(LAPACKE::${_comp} ${LIB_TYPE} IMPORTED GLOBAL)
      if(LAPACKE_INCLUDE_DIRS)
        set_target_properties(LAPACKE::${_comp} PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS}")
      endif()
      if(EXISTS "${LAPACKE_${_comp}_LIB}" AND NOT "${LIB_TYPE}" STREQUAL INTERFACE)
          set_target_properties(LAPACKE::${_comp} PROPERTIES
            IMPORTED_LOCATION "${LAPACKE_${_comp}_LIB}")
        if(CMAKE_Fortran_COMPILER_WORKS)
          set_target_properties(LAPACKE::${_comp} PROPERTIES
            LINKER_LANGUAGE "Fortran")
        endif()
      endif()
      list(APPEND _fallback_${_comp} "${MATH_LIB}")
      set_target_properties(LAPACKE::${_comp} PROPERTIES
          INTERFACE_LINK_LIBRARIES "${_fallback_${_comp}}")
    endif()
  endforeach()
  target_link_libraries(LAPACKE::LAPACKE INTERFACE LAPACKE::LAPACK LAPACKE::CBLAS LAPACKE::BLAS)

  # ----------------------------------------------------------------------------

  if(NOT LAPACKE_FIND_QUIETLY)
    message(STATUS "Found LAPACKE and defined the following imported targets:")
    foreach(_comp ${LAPACKE_FIND_COMPONENTS})
      message(STATUS "  - LAPACKE::${_comp}:")
      message(STATUS "      + include:      ${LAPACKE_INCLUDE_DIRS}")
      message(STATUS "      + library:      ${LAPACKE_${_comp}_LIB}")
      message(STATUS "      + dependencies: ${_fallback_${_comp}}")
    endforeach()
  endif()
endif()

# ==============================================================================

mark_as_advanced(
  LAPACKE_FOUND
  LAPACKE_INCLUDE_DIRS
  LAPACKE_LIBRARIES
  )
