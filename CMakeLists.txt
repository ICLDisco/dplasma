cmake_minimum_required (VERSION 3.15)
project (DPLASMA C)

include(CMakeDependentOption)
include(CMakePushCheckState)
include(CheckSymbolExists)
include(CheckCSourceCompiles)
include(CheckIncludeFiles)
include(GNUInstallDirs)

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
if(POLICY CMP0079)
  # CMP0079: CMake 3.13 now allows the target_link_libraries() command to be called from any directory
  # to add link dependencies and link interface libraries to targets created in other directories.
  cmake_policy(SET CMP0079 NEW)
endif(POLICY CMP0079)
if(POLICY CMP0094)
  # CMP0094: Starting with CMake 3.16, all FIND_Python will use the first matching version instead of
  #          of searching for the largest available version number (which defeats our selection logic)
  cmake_policy(SET CMP0094 NEW)
ENDIF(POLICY CMP0094)


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

option(BUILD_64bits
  "Build 64 bits mode" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_SHARED_LIBS
       "Enable building shared DPLASMA lib (default ON)" ON)

option(SUPPORT_FORTRAN
       "Enable support for Fortran bindings (default ON)" ON)
# In CMake 3.12, enable_language OPTIONAL has no effect and will still
# enable the language for a compiler that does not work/is not present
# causing all sort of problems downstream. For now use CheckLanguage before.
include(CheckLanguage)
if( SUPPORT_FORTRAN )
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
  endif()
endif( SUPPORT_FORTRAN )

option(SUPPORT_C11
       "Enable support for C11 capabilities. Might not work with OpenMP support due to an OpenMP compilers restriction (default ON)." ON)
if( SUPPORT_C11 )
  set(CMAKE_C_STANDARD 11)
else (SUPPORT_C11)
  set(CMAKE_C_STANDARD 99)
endif( SUPPORT_C11 )
set(CMAKE_C_STANDARD_REQUIRED ON)

option(DPLASMA_WITH_SCALAPACK_WRAPPER
  "Build a ScaLAPACK compatible API to DPLASMA" ON)

option(DPLASMA_WITH_RECURSIVE
  "Enable recursive kernels to be called when available" ON)

option(DPLASMA_TRACE_KERNELS
  "Enable tracing dplasma kernel calls for debugging (slow)" OFF)

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
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

#
# Find packages to decide if we compile DPLASMA support or not
#

include(AddDocumentedFiles)

if (CMAKE_Fortran_COMPILER_WORKS)
include(FortranCInterface)
FortranCInterface_HEADER(src/include/dplasma_f2c.h
                         MACRO_NAMESPACE "DPLASMA_F2C_"
                         SYMBOL_NAMESPACE "dplasma_f2c_")
else()
# No Fortran compiler to test mangling. So assume trailing_ and hope for the best.
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/src/include/dplasma_f2c.h
"
#ifndef DPLASMA_F2C_HEADER_INCLUDED
#define DPLASMA_F2C_HEADER_INCLUDED
#define DPLASMA_F2C_GLOBAL(name,NAME) name##_
#endif
")
endif()

#set(BLA_STATIC OFF)

set(BLA_VENDOR CACHE STRING "Possible valid values are: Intel10_64lp_seq, IBMESSL, Apple, or All, which will default to something like Goto or ATLAS")
# Goto, ATLAS, PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_64lp_seq, ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic
set_property(CACHE BLA_VENDOR PROPERTY STRINGS Intel10_64lp_seq IBMESSL Apple Goto ATLAS)
find_package(LAPACKE REQUIRED)

check_include_files(complex.h DPLASMA_HAVE_COMPLEX_H)

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
if( APPLE )
  CMAKE_PUSH_CHECK_STATE()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES}")
  include(CheckFunctionExists)
  check_function_exists(appleblas_sgeadd PARSEC_OSX_VECLIB)
  CMAKE_POP_CHECK_STATE()
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
      MESSAGE(WARNING "Detected Apple vecLIB! the Accelerate Framework has a non-standard BLAS interface that causes DPLASMA to produce wrong results in single precision. Restricting DPLASMA to ${NEW_PREC} (from ${DPLASMA_PRECISIONS})")
      set(DPLASMA_PRECISIONS ${NEW_PREC})
    endif( prec_change )
  endif( PARSEC_OSX_VECLIB )
endif( APPLE )

#
# Find or Build PaRSEC package
option(DPLASMA_FORCE_INTERNAL_PARSEC "Force using the internal PaRSEC submodule" OFF)
if(NOT DPLASMA_FORCE_INTERNAL_PARSEC AND EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
    # We were not instructed to ignore external PaRSEC during the fist configuration, but we did
    # use the submodule in the past. Keep doing so now, except if a PaRSEC_ROOT has been supplied
    # explicit, in which case we will be noisy about this being unusual.
    if(PaRSEC_ROOT)
      message(WARNING "This source directory contains the PaRSEC submodule, but PaRSEC_ROOT=${PaRSEC_ROOT} has been supplied.\n  This build will now ignore the submodule and use the external PaRSEC in ${PaRSEC_ROOT}.")
    else()
      # PaRSEC set to autodetect, and we used the submodule in the past, keep doing so.
      # Not in the cache, we want to do the check for PaRSEC_ROOT again next time.
      set(DPLASMA_FORCE_INTERNAL_PARSEC ON)
    endif()
endif()
# Find external PaRSEC if not disallowed by submodule
if(NOT DPLASMA_FORCE_INTERNAL_PARSEC)
    find_package(PaRSEC HINTS "${PaRSEC_ROOT}" QUIET)
endif()

if(NOT TARGET PaRSEC::parsec AND NOT TARGET PaRSEC::parsec_ptgpp)
  if(PaRSEC_ROOT AND DPLASMA_FORCE_INTERNAL_PARSEC)
      message(FATAL_ERROR "User requested PaRSEC with PaRSEC_ROOT=${PaRSEC_ROOT} and DPLASMA_FORCE_INTERNAL_PARSEC is set; this is incoherent")
  elseif(PaRSEC_ROOT)
    message(FATAL_ERROR "User requested PaRSEC with PaRSEC_ROOT=${PaRSEC_ROOT} not found")
  endif()
  # we have to build it
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
      # Update submodules as needed
      option(GIT_SUBMODULE "Check submodules during build" ON)
      if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --force --init --recursive
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
    # Import targets and load dependent libraries
    find_package(parsec-submodule NO_DEFAULT_PATH PATHS ${CMAKE_CURRENT_BINARY_DIR}/parsec REQUIRED)
    set(DPLASMA_FORCE_INTERNAL_PARSEC ${PaRSEC_FOUND} CACHE BOOL "Force using the internal PaRSEC submodule (submodule used in previous configure)" FORCE)
    message(STATUS "#             Configuring internal submodule PaRSEC runtime: DONE!")
    message(STATUS "########################################################################")
  endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/parsec/CMakeLists.txt")
endif(NOT TARGET PaRSEC::parsec AND NOT TARGET PaRSEC::parsec_ptgpp)


option(DPLASMA_HAVE_CUDA "Use CUDA to accelerate DPLASMA routines" ${PARSEC_HAVE_CUDA})
if(DPLASMA_HAVE_CUDA)
  if(NOT TARGET CUDA::cusolver)
    find_package(CUDAToolkit REQUIRED)
  endif(NOT TARGET CUDA::cusolver)
  message(STATUS "CUDA support for DPLASMA enabled")
endif()


set(CMAKE_Fortran_LINKER_PREFERENCE_PROPAGATES TRUE)
add_subdirectory(src)
add_subdirectory(src/scalapack_wrappers)
IF( BUILD_TESTING )
  add_subdirectory(tests)
ENDIF( BUILD_TESTING )
add_subdirectory(examples)
add_subdirectory(docs)


set(DPLASMA_SHARE_DIRECTORY_SOURCE ${PROJECT_SOURCE_DIR}/share)
set(DPLASMA_SHARE_DIRECTORY_INSTALL ${CMAKE_INSTALL_DATADIR}/dplasma)

# Configuration header
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/src/include/dplasma/config.h.in"
  "${PROJECT_BINARY_DIR}/src/include/dplasma/config.h")
install(FILES
  "${PROJECT_BINARY_DIR}/src/include/dplasma/config.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dplasma)

# Export Targets
install(EXPORT dplasma-targets # rules to export the targets built here
  FILE dplasma-targets.cmake
  NAMESPACE DPLASMA::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/dplasma)

install(FILES share/help-dplasma.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/dplasma)

# cmake config files
include(CMakePackageConfigHelpers)

configure_package_config_file(cmake_modules/dplasma-config.cmake.in
  ${PROJECT_BINARY_DIR}/cmake_modules/dplasma-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/dplasma
  PATH_VARS
    CMAKE_INSTALL_INCLUDEDIR
    CMAKE_INSTALL_DATADIR
    CMAKE_INSTALL_LIBDIR
    CMAKE_INSTALL_BINDIR
    CMAKE_INSTALL_LIBEXECDIR)

write_basic_package_version_file(${PROJECT_BINARY_DIR}/cmake_modules/dplasma-config-version.cmake
    VERSION ${DPLASMA_VERSION_MAJOR}.${DPLASMA_VERSION_MINOR}
    COMPATIBILITY SameMajorVersion)

install(FILES
  ${PROJECT_BINARY_DIR}/cmake_modules/dplasma-config.cmake
  ${PROJECT_BINARY_DIR}/cmake_modules/dplasma-config-version.cmake
  ${PROJECT_SOURCE_DIR}/cmake_modules/FindLAPACKE.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/dplasma)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/src/include/dplasma.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dplasma.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dplasma.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

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
