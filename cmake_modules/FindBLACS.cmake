#.rst:
# FindBLACS
# -------------
#
# Find BLACS include dirs and libraries
#
# Use this module by invoking find_package with the form::
#
#   find_package(BLACS
#     [REQUIRED]             # Fail with error if BLACS is not found
#     )
#
# This module defines::
#
#   BLACS_FOUND            - True if headers and requested libraries were found
#   BLACS_INCLUDE_DIRS     - BLACS include directories
#   BLACS_LIBRARIES        - BLACS component libraries to be linked
#
#
# This module reads hints about search locations from variables
# (either CMake variables or environment variables)::
#
#   BLACS_ROOT             - Preferred installation prefix for BLACS
#   BLACS_DIR              - Preferred installation prefix for BLACS
#
#
# The following :prop_tgt:`IMPORTED` targets are also defined::
#
#   BLACS::blacs       - Imported target for the BLACS library
#
# ==============================================================================
# Copyright (c) 2023      The University of Tennessee and The University
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

if(NOT TARGET LAPACKE::LAPACK OR NOT TARGET LAPACKE::BLAS)
  find_package(LAPACKE COMPONENTS LAPACK;BLAS QUIET)
endif()
if(NOT TARGET MPI::MPI_C)
  find_package(MPI QUIET)
endif()

# Is BLACS integrated within BLAS/LAPACK?
cmake_push_check_state()
set(CMAKE_REQUIRED_LIBRARIES LAPACKE::LAPACK;LAPACKE::BLAS;MPI::MPI_C)
include(CheckCSourceCompiles)
check_c_source_compiles("extern void Cblacs2sys_handle(void); int main(void) { Cblacs2sys_handle(); return 0; }" BLAS_HAS_2SYSHANDLE)
cmake_pop_check_state()

if(BLAS_HAS_BLACS)
list(APPEND BLACS_REQUIRED_VARS BLAS_HAS_2SYSHANDLE)

else(BLAS_HAS_BLACS)
  if(BLAS_LIBRARIES)
    list(GET BLAS_LIBRARIES 0 blas_main_library)
    get_filename_component(blas_dir ${blas_main_library} DIRECTORY)
  endif()
  set(BLACS_SEARCH_PATHS
    ${blas_dir}
    $ENV{MKLROOT}
    $ENV{MKLROOT}/mkl
    /usr/local/opt # homebrew on mac
    /opt
    /opt/local
  )

  set(LIB_PATH_SUFFIXES
    intel64
    mkl/lib mkl/lib/intel64
    compiler/lib compiler/lib/intel64
  )

  if(APPLE)
    list(APPEND LIB_PATH_SUFFIXES scalapack/lib openblas/lib)
  elseif(WIN32)
    list(APPEND BLACS_SEARCH_PATHS "C:/Program Files (x86)/ScaLAPACK")
    list(APPEND BLACS_SEARCH_PATHS "C:/Program Files/ScaLAPACK")
  endif()

  # We actually need only the BLACS part, but it is inside ref-scalapack
  # We use the lp64 always, as this is how our own generic headers declare the functions
  find_library(BLACS_LIBRARY
    NAMES mkl_blacs_openmpi_lp64 mkl_blacs_intelmpi_lp64 scalapack
    NAMES_PER_DIR
    PATHS ${BLACS_SEARCH_PATHS}
    PATH_SUFFIXES ${LIB_PATH_SUFFIXES})
  list(APPEND BLACS_REQUIRED_VARS BLACS_LIBRARY)

  if(BLACS_LIBRARY)
    # check it works
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES ${BLACS_LIBRARY};LAPACKE::LAPACK;LAPACKE::BLAS;MPI::MPI_C)
    include(CheckCSourceCompiles)
    check_c_source_compiles("extern void Cblacs2sys_handle(void); int main(void) { Cblacs2sys_handle(); return 0; }" BLACS_HAS_2SYSHANDLE)
    cmake_pop_check_state()
    list(APPEND BLACS_REQUIRED_VARS BLACS_HAS_2SYSHANDLE)
  endif()
endif(BLAS_HAS_BLACS)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLACS
  FOUND_VAR BLACS_FOUND
  REQUIRED_VARS ${BLACS_REQUIRED_VARS})

if(BLACS_FOUND)
  set(BLACS_INCLUDE_DIRS "")
  set(BLACS_LIBRARIES "${BLACS_LIBRARY}")
  get_filename_component(LIB_EXT "${BLACS_LIBRARY}" EXT)
  if(LIB_EXT STREQUAL "" OR LIB_EXT STREQUAL ".framework")
     set(LIB_TYPE INTERFACE)
  elseif(LIB_EXT STREQUAL ".a" OR LIB_EXT STREQUAL ".lib")
     set(LIB_TYPE STATIC)
  else()
    set(LIB_TYPE SHARED)
  endif()
  add_library(BLACS::blacs ${LIB_TYPE} IMPORTED GLOBAL)
  if(BLACS_INCLUDE_DIRS)
    set_target_properties(BLACS::blacs PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${BLACS_INCLUDE_DIRS}")
  endif()
  if(EXISTS "${BLACS_LIBRARY}" AND NOT "${LIB_TYPE}" STREQUAL INTERFACE)
    set_target_properties(BLACS::blacs PROPERTIES
      IMPORTED_LOCATION "${BLACS_LIBRARY}")
  endif()
  target_link_libraries(BLACS::blacs INTERFACE LAPACKE::LAPACK LAPACKE::BLAS)
  mark_as_advanced(BLACS_INCLUDE_DIRS BLACS_LIBRARIES)
endif()
