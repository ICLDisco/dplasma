#.rst:
# FindScaLAPACK
# -------------
#
# Find ScaLAPACK include dirs and libraries
#
# Use this module by invoking find_package with the form::
#
#   find_package(ScaLAPACK
#     [REQUIRED]             # Fail with error if ScaLAPACK is not found
#     )
#
# This module defines::
#
#   ScaLAPACK_FOUND            - True if headers and requested libraries were found
#   ScaLAPACK_INCLUDE_DIRS     - ScaLAPACK include directories
#   ScaLAPACK_LIBRARIES        - ScaLAPACK component libraries to be linked
#
#
# This module reads hints about search locations from variables
# (either CMake variables or environment variables)::
#
#   ScaLAPACK_ROOT             - Preferred installation prefix for ScaLAPACK
#   ScaLAPACK_DIR              - Preferred installation prefix for ScaLAPACK
#
#
# The following :prop_tgt:`IMPORTED` targets are also defined::
#
#   ScaLAPACK::scalapack       - Imported target for the ScaLAPACK library
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

# Is ScaLAPACK integrated within BLAS/LAPACK?
cmake_push_check_state()
set(CMAKE_REQUIRED_LIBRARIES LAPACKE::LAPACK;LAPACKE::BLAS;MPI::MPI_C)
include(CheckCSourceCompiles)
check_c_source_compiles("extern void Cblacs2sys_handle(void); int main(void) { Cblacs2sys_handle(); return 0; }" BLAS_HAS_ScaLAPACK)
cmake_pop_check_state()

if(BLAS_HAS_ScaLAPACK)
  list(APPEND ScaLAPACK_REQUIRED_VARS BLAS_HAS_ScaLAPACK)

else(BLAS_HAS_ScaLAPACK)
  if(BLAS_LIBRARIES)
    list(GET BLAS_LIBRARIES 0 blas_main_library)
    get_filename_component(blas_dir ${blas_main_library} DIRECTORY)
  endif()
  set(ScaLAPACK_SEARCH_PATHS
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
    list(APPEND ScaLAPACK_SEARCH_PATHS "C:/Program Files (x86)/ScaLAPACK")
    list(APPEND ScaLAPACK_SEARCH_PATHS "C:/Program Files/ScaLAPACK")
  endif()

  # We actually need only the BLACS part, but it is inside ref-scalapack
  # We use the lp64 always, as this is how our own generic headers declare the functions
  find_library(ScaLAPACK_LIBRARY
    NAMES mkl_blacs_openmpi_lp64 mkl_blacs_intelmpi_lp64 scalapack
    NAMES_PER_DIR
    PATHS ${ScaLAPACK_SEARCH_PATHS}
    PATH_SUFFIXES ${LIB_PATH_SUFFIXES})
  list(APPEND ScaLAPACK_REQUIRED_VARS ScaLAPACK_LIBRARY)

  if(ScaLAPACK_LIBRARY)
    # check it works
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES ${ScaLAPACK_LIBRARY};LAPACKE::LAPACK;LAPACKE::BLAS;MPI::MPI_C)
    include(CheckCSourceCompiles)
    check_c_source_compiles("extern void Cblacs2sys_handle(void); int main(void) { Cblacs2sys_handle(); return 0; }" ScaLAPACK_FOUND_BLACS)
    cmake_pop_check_state()
    list(APPEND ScaLAPACK_REQUIRED_VARS ScaLAPACK_FOUND_BLACS)
  endif()
endif(BLAS_HAS_ScaLAPACK)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScaLAPACK
  FOUND_VAR ScaLAPACK_FOUND
  REQUIRED_VARS ${ScaLAPACK_REQUIRED_VARS})

if(ScaLAPACK_FOUND)
  set(ScaLAPACK_INCLUDE_DIRS "")
  set(ScaLAPACK_LIBRARIES "${ScaLAPACK_LIBRARY}")
  get_filename_component(LIB_EXT "${ScaLAPACK_LIBRARY}" EXT)
  if(LIB_EXT STREQUAL "" OR LIB_EXT STREQUAL ".framework")
     set(LIB_TYPE INTERFACE)
  elseif(LIB_EXT STREQUAL ".a" OR LIB_EXT STREQUAL ".lib")
     set(LIB_TYPE STATIC)
  else()
    set(LIB_TYPE SHARED)
  endif()
  add_library(ScaLAPACK::scalapack ${LIB_TYPE} IMPORTED GLOBAL)
  if(ScaLAPACK_INCLUDE_DIRS)
    set_target_properties(ScaLAPACK::scalapack PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${ScaLAPACK_INCLUDE_DIRS}")
  endif()
  if(EXISTS "${ScaLAPACK_LIBRARY}" AND NOT "${LIB_TYPE}" STREQUAL INTERFACE)
    set_target_properties(ScaLAPACK::scalapack PROPERTIES
      IMPORTED_LOCATION "${ScaLAPACK_LIBRARY}")
    if(CMAKE_Fortran_COMPILER_WORKS)
      set_target_properties(ScaLAPACK::scalapack PROPERTIES
         LINKER_LANGUAGE "Fortran")
    endif()
  endif()
  target_link_libraries(ScaLAPACK::scalapack INTERFACE LAPACKE::LAPACK LAPACKE::BLAS)
  mark_as_advanced(ScaLAPACK_INCLUDE_DIRS ScaLAPACK_LIBRARIES)
endif()
