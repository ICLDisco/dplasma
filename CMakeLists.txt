cmake_minimum_required (VERSION 3.12)
project (DPLASMA C)

include(CMakeDependentOption)
include(CMakePushCheckState)
include(CheckSymbolExists)
include(CheckCSourceCompiles)

# The current version number
set (DPLASMA_VERSION_MAJOR 2)
set (DPLASMA_VERSION_MINOR 0)
set (DPLASMA_VERSION_PATCH 0)

SET(DPLASMA_VERSION "${DPLASMA_VERSION_MAJOR}.${DPLASMA_VERSION_MINOR}")

set(CMAKE_NO_SYSTEM_FROM_IMPORTED True)

# CMake Policies Tuning
if(POLICY CMP0074)
  # CMP0074: Starting with CMake 3.12, all FIND_<something> use <something>_ROOT in the search path
  #          in addition to the specified paths
  cmake_policy(SET CMP0074 NEW)
ENDIF(POLICY CMP0074)

# CTest system
SET(DART_TESTING_TIMEOUT 120)
enable_testing()
include(CTest)

set(CTEST_SHM_LAUNCHER
      "" CACHE STRING
      "A command to run shared memory testings")
string(REPLACE " " ";" SHM_TEST_CMD_LIST "${CTEST_SHM_LAUNCHER}")
set(CTEST_MPI_LAUNCHER
      "mpiexec -n" CACHE STRING
      "A command to run distributed memory testings")
if( "${CTEST_MPI_LAUNCHER}" STREQUAL "" )
    MESSAGE(WARNING "MPI tests will most likely not work: 'CTEST_MPI_LAUNCHER' is not set")
endif()
string(REPLACE " " ";" MPI_TEST_CMD_LIST "${CTEST_MPI_LAUNCHER}")
set(CTEST_GPU_LAUNCHER_OPTIONS
      "" CACHE STRING
      "Options to pass to the CTEST_MPI_LAUNCHER to select GPU nodes")

# DPLASMA specific option.
option(BUILD_SHARED_LIBS
       "Enable building shared DPLASMA lib (default ON)" ON)

option(SUPPORT_FORTRAN
       "Enable support for Fortran bindings (default ON)" ON)
if( SUPPORT_FORTRAN )
  enable_language(Fortran OPTIONAL)
endif( SUPPORT_FORTRAN )

option(DPLASMA_WITH_RECURSIVE
  "Enable recursive kernels to be called when available" OFF)

mark_as_advanced( DPLASMA_GPU_WITH_MAGMA )
option(DPLASMA_GPU_WITH_MAGMA
  "Enable GPU support using MAGMA kernels" OFF)
if(DPLASMA_GPU_WITH_MAGMA)
  message(WARNING "MAGMA is not supported yet, ignored.")
endif()

mark_as_advanced(DPLASMA_DEBUG_QR_PIVGEN)
option(DPLASMA_DEBUG_QR_PIVGEN
  "Enable the QR pivgen testings" OFF)

STRING(COMPARE EQUAL ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} DPLASMA_BUILD_INPLACE)
#
# Find packages to decide if we compile DPLASMA support or not
#
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

include(AddDocumentedFiles)

include(FortranCInterface)
FortranCInterface_HEADER(include/dplasma_f2c.h
                         MACRO_NAMESPACE "DPLASMA_F2C_"
                         SYMBOL_NAMESPACE "dplasma_f2c_")

#set(BLA_STATIC OFF)

set(BLA_VENDOR CACHE STRING "Possible valid values are: Intel10_64lp_seq, IBMESSL, Apple, or All, which will default to something like Goto or ATLAS")
# Goto, ATLAS, PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_64lp_seq, ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic
set_property(CACHE BLA_VENDOR PROPERTY STRINGS Intel10_64lp_seq IBMESSL Apple Goto ATLAS)
find_package(LAPACKE REQUIRED)

check_symbol_exists(timersub "sys/time.h" DPLASMA_TIMERSUB_NEED__GNU_SOURCE)
if( NOT DPLASMA_TIMERSUB_NEED__GNU_SOURCE )
  unset(DPLASMA_TIMERSUB_NEED__GNU_SOURCE CACHE)
  CMAKE_PUSH_CHECK_STATE()
  list(APPEND CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE")
  check_symbol_exists(timersub "sys/time.h" DPLASMA_TIMERSUB_NEED__GNU_SOURCE)
  CMAKE_POP_CHECK_STATE()
  if( DPLASMA_TIMERSUB_NEED__GNU_SOURCE )
    add_definitions(-D_GNU_SOURCE)
  endif( DPLASMA_TIMERSUB_NEED__GNU_SOURCE )
endif( NOT DPLASMA_TIMERSUB_NEED__GNU_SOURCE )
set(DPLASMA_NEED__GNU_SOURCE ${DPLASMA_TIMERSUB_NEED__GNU_SOURCE} CACHE BOOL "True if DPLASMA needs _GNU_SOURCE")

check_symbol_exists(asprintf "stdio.h" DPLASMA_ASPRINTF_NEED__GNU_SOURCE)
if( NOT DPLASMA_ASPRINTF_NEED__GNU_SOURCE )
  unset(DPLASMA_ASPRINTF_NEED__GNU_SOURCE CACHE)
  CMAKE_PUSH_CHECK_STATE()
  list(APPEND CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE")
  check_symbol_exists(asprintf "stdio.h" DPLASMA_ASPRINTF_NEED__GNU_SOURCE)
  CMAKE_POP_CHECK_STATE()
  if( DPLASMA_ASPRINTF_NEED__GNU_SOURCE )
    set(DPLASMA_NEED__GNU_SOURCE True CACHE BOOL "True if DPLASMA needs _GNU_SOURCE")
    add_definitions(-D_GNU_SOURCE)
  endif( DPLASMA_ASPRINTF_NEED__GNU_SOURCE )
endif( NOT DPLASMA_ASPRINTF_NEED__GNU_SOURCE )

if( DPLASMA_GPU_WITH_MAGMA )
  find_package(MAGMA)
  if( MAGMA_FOUND )
  else( MAGMA_FOUND )
  endif( MAGMA_FOUND )
endif( DPLASMA_GPU_WITH_MAGMA )

if(NOT DPLASMA_PRECISIONS)
  set(DPLASMA_PRECISIONS "s;d;c;z" CACHE STRING "The precisions to compile in dplasma (accepts a colon separated list of s;d;c;z)" FORCE)
else()
  set(DPLASMA_PRECISIONS "${DPLASMA_PRECISIONS}" CACHE STRING "The precisions to compile in dplasma (accepts a colon separated list of s;d;c;z)" FORCE)
endif()
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tools/PrecisionGenerator/")
Add_Subdirectory(tools/PrecisionGenerator)

# Veclib on OSX diverges from the standard BLAS interface. As a result, single
# precision dplasma tests (real and complex) deliver wrong results.
# Until we have a fix, make sure to disable these 2 presicisions with vecLib.
if( PARSEC_OSX )
  check_function_exists(appleblas_sgeadd PARSEC_OSX_VECLIB)
  if( PARSEC_OSX_VECLIB )
    set(tmp_fixed_prec "")
    foreach(prec ${DPLASMA_PRECISIONS})
      STRING(COMPARE EQUAL "${prec}" "s" is_sprec)
      STRING(COMPARE EQUAL "${prec}" "c" is_cprec)
      if ( NOT (is_sprec OR is_cprec) )
        set(tmp_fixed_prec "${tmp_fixed_prec} ${prec}")
      endif ( NOT (is_sprec OR is_cprec) )
    endforeach()
    STRING(STRIP "${tmp_fixed_prec}" tmp_fixed_prec)
    STRING(REPLACE " " ";" NEW_PREC "${tmp_fixed_prec}")
    STRING(COMPARE NOTEQUAL "${NEW_PREC}" "${DPLASMA_PRECISIONS}" prec_change)
    if( prec_change )
        MESSAGE(WARNING "Detected Apple vecLIB! the Accelerate Framework has a non-standard BLAS interface that causes Dplasma to produce wrong results in single precision. Restricting DPLASMA to ${NEW_PREC} (from ${DPLASMA_PRECISIONS})")
      set(DPLASMA_PRECISIONS ${NEW_PREC})
    endif( prec_change )
  endif( PARSEC_OSX_VECLIB )
endif( PARSEC_OSX )

set(DPLASMA_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/dplasma"
                           CACHE STRING "DPLASMA Installation Prefix")
#
# Find or Build PaRSEC package
STRING(COMPARE EQUAL "-internal-" "${PaRSEC_ROOT}" DPLASMA_FORCE_INTERNAL)
STRING(COMPARE NOTEQUAL "" "${PaRSEC_ROOT}" PARSEC_EXTERNAL_REQUIRED)
if(NOT DPLASMA_FORCE_INTERNAL)
  find_package(PaRSEC HINTS "${PaRSEC_ROOT}" PATH_SUFFIXES share/cmake/parsec cmake_modules QUIET)
endif()

if(NOT TARGET PaRSEC::parsec AND NOT TARGET PaRSEC::parsec_ptgpp)
  if(PARSEC_EXTERNAL_REQUIRED)
    message(FATAL_ERROR "User requested PaRSEC not found with PaRSEC_ROOT=${PaRSEC_ROOT}")
  endif()
  # we have to build it
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
      # Update submodules as needed
      option(GIT_SUBMODULE "Check submodules during build" ON)
      if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --force --init --recursive --remote
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
          message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
      endif()
    endif()
  endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")

  # TODO: submodule stuff from dplasma
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
    message(FATAL_ERROR "The PaRSEC submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
  else(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
    # options for PaRSEC?
    message(STATUS "########################################################################")
    message(STATUS "#             Configuring internal submodule PaRSEC runtime!")
    # Thanks to CMake policy CMP0077 we cannot set a variable before it is used with option(). Thus,
    # we need to make it an option here, so that cmake accept it and does the right thing (use the value
    # we need and not the default value provided to the first OPTION call.
    option(PARSEC_WITH_DEVEL_HEADERS "" ON)
    add_subdirectory(parsec)
    # that's the limit of mixing namespace, export, and building a project, you're stuck with mismatching targets
    add_library(PaRSEC::parsec ALIAS parsec)
    add_executable(PaRSEC::parsec_ptgpp ALIAS parsec_ptgpp)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/parsec/cmake_modules)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/parsec
                        ${CMAKE_CURRENT_SOURCE_DIR}/parsec/parsec/include)
    include_directories($<$<NOT:${DPLASMA_BUILD_INPLACE}>:${CMAKE_CURRENT_BINARY_DIR}/parsec
                                                          ${CMAKE_CURRENT_BINARY_DIR}/parsec/parsec/include>)
    # I don't understand why CMake does not propagate the parsec target include requirements
    # but it does not. So we need to forcefully add the MPI C header path.
    include_directories(${MPI_C_HEADER_DIR})
    message(STATUS "#             Configuring internal submodule PaRSEC runtime: DONE!")
    message(STATUS "########################################################################")
  endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
endif(NOT TARGET PaRSEC::parsec AND NOT TARGET PaRSEC::parsec_ptgpp)


if(PARSEC_HAVE_CUDA)
  # CUDA is pulled by PaRSEC
  set(DPLASMA_HAVE_CUDA ${PARSEC_HAVE_CUDA} CACHE BOOL "DPLASMA compiled with CUDA support")
  message(STATUS "CUDA support for DPLASMA enabled")
endif(PARSEC_HAVE_CUDA)

include (CheckIncludeFiles)
check_include_files(complex.h DPLASMA_HAVE_COMPLEX_H)

if(BLAS_FOUND)
  set(CMAKE_Fortran_LINKER_PREFERENCE_PROPAGATES TRUE)
  Add_Subdirectory(include)
  Add_Subdirectory(cores)
  Add_Subdirectory(lib)
  IF( BUILD_TESTING )
    Add_Subdirectory(tests)
  ENDIF( BUILD_TESTING )
  Add_Subdirectory(examples)
endif(BLAS_FOUND)

Add_Subdirectory(docs)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/dplasma.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dplasma.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dplasma.pc"
        DESTINATION ${DPLASMA_INSTALL_PREFIX}/lib/pkgconfig)

# Prepare the contrib/build_with_dplasma directory
configure_file (
  contrib/build_with_dplasma/Makefile.in
  ${CMAKE_CURRENT_BINARY_DIR}/contrib/build_with_dplasma/Makefile @ONLY)

file(GLOB_RECURSE pattern_files RELATIVE
     "${CMAKE_CURRENT_SOURCE_DIR}/" "contrib/build_with_dplasma/*.[c|h]"
                                    "contrib/build_with_dplasma/*.jdf")
foreach( pattern_file ${pattern_files} )
  if(NOT DPLASMA_BUILD_INPLACE)
    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${pattern_file}"
      COMMAND ${CMAKE_COMMAND} -E copy
              "${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file}"
              "${CMAKE_CURRENT_BINARY_DIR}/${pattern_file}"
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file}"
    )
  endif(NOT DPLASMA_BUILD_INPLACE)
  list( APPEND pattern_files_dest "${pattern_file}" )
endforeach( pattern_file )
add_custom_target(build_with_dplasma ALL DEPENDS ${pattern_files_dest})
