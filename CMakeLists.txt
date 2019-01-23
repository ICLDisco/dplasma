cmake_minimum_required (VERSION 3.0)
project (DPLASMA C)

include(CMakeDependentOption)
include(CMakePushCheckState)

# The current version number
set (DPLASMA_VERSION_MAJOR 2)
set (DPLASMA_VERSION_MINOR 0)
set (DPLASMA_VERSION_PATCH 0)

SET(DPLASMA_VERSION "${PARSEC_VERSION_MAJOR}.${PARSEC_VERSION_MINOR}")

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
option(DPLASMA_WITH_RECURSIVE
  "Enable recursive kernels to be called when available" OFF)
if(DPLASMA_WITH_RECURSIVE)
  # TODO: This will need to be added to the futur dplasma_config.h.in
  add_definitions(-DPARSEC_HAVE_RECURSIVE)
endif()

mark_as_advanced( PARSEC_GPU_WITH_MAGMA )
option(DPLASMA_GPU_WITH_MAGMA
  "Enable GPU support using MAGMA kernels" OFF)
if(DPLASMA_GPU_WITH_MAGMA)
  message(WARNING "MAGMA is not supported yet, ignored.")
endif()

mark_as_advanced(DPLASMA_DEBUG_QR_PIVGEN)
option(DPLASMA_DEBUG_QR_PIVGEN
  "Enable the QR pivgen testings" OFF)

#
# Find or Build PaRSEC package
#
list(APPEND CMAKE_MODULE_PATH /home/dgenet/pr/split/install/parsec-alone/lib/cmake)

#list(APPEND PKG_CONFIG_PATH /home/dgenet/pr/split/install/parsec-alone/lib/pkgconfig)
message(STATUS "XXX1/ ${CMAKE_MODULE_PATH}")

find_package(PaRSEC)

if (PARSEC_ptgpp_EXE AND PARSEC_LIBRARIES)
    # continue and hope for the best
message(STATUS "I found that ${PARSEC_ptgpp_EXE} and that ${PARSEC_LIBRARIES}")
else (PARSEC_ptgpp_EXE AND PARSEC_LIBRARIES)
    message(STATUS "I didnt find PARSEC, I will have to build it")
    # we have to build whatever is in parsec/
    include(parsec/CMakeLists.txt)
endif (PARSEC_ptgpp_EXE AND PARSEC_LIBRARIES)

#
# Find packages to decide if we compile DPLASMA support or not
#
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

include(FortranCInterface)

find_package(COREBLAS)
if( PLASMA_SRC )
   list(APPEND COREBLAS_INCLUDE_DIRS ${PLASMA_SRC})
endif( PLASMA_SRC )

if( DPLASMA_GPU_WITH_MAGMA )
  find_package(MAGMA)
  if( MAGMA_FOUND )
  else( MAGMA_FOUND )
  endif( MAGMA_FOUND )
endif( DPLASMA_GPU_WITH_MAGMA )

include_directories(${COREBLAS_INCLUDE_DIRS})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if( ${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
else( ${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/include)
endif( ${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )

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

if(COREBLAS_FOUND)
  Add_Subdirectory(include)
  Add_Subdirectory(cores)
  Add_Subdirectory(lib)
  IF( BUILD_TESTING )
    Add_Subdirectory(testing)
  ENDIF( BUILD_TESTING )
  Add_Subdirectory(examples)
endif(COREBLAS_FOUND)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/dplasma.pc.in"
  "${PROJECT_INCLUDE_DIR}/dplasma.pc" @ONLY)
install(FILES "${PROJECT_INCLUDE_DIR}/dplasma.pc"
        DESTINATION ${DPLASMA_INSTALL_PREFIX}/lib/pkgconfig)

# Prepare the contrib/build_with_parsec directory
configure_file (
  contrib/build_with_parsec/Makefile.in
  ${CMAKE_CURRENT_BINARY_DIR}/contrib/build_with_parsec/Makefile @ONLY)

file(GLOB_RECURSE pattern_files RELATIVE
     "${CMAKE_CURRENT_SOURCE_DIR}/" "contrib/build_with_parsec/*.[c|h]"
                                    "contrib/build_with_parsec/*.jdf")
foreach( pattern_file ${pattern_files} )
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${pattern_file}"
    COMMAND cmake -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file}"
            "${CMAKE_CURRENT_BINARY_DIR}/${pattern_file}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file}"
  )
  list( APPEND pattern_files_dest "${pattern_file}" )
endforeach( pattern_file )
add_custom_target(build_with_parsec ALL DEPENDS ${pattern_files_dest})
