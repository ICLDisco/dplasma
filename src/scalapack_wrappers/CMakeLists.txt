# Generate ScaLAPACK wrappers provided by DPLASMA
include(DplasmaCompilerFlags)

if(NOT DPLASMA_WITH_SCALAPACK_WRAPPER)
  return() # told to not build them
endif(NOT DPLASMA_WITH_SCALAPACK_WRAPPER)

if(NOT TARGET BLACS::blacs)
  message(WARNING "DPLASMA ScaLAPACK wrappers were requested but BLACS could not be found: disabling the wrappers")
  return()
endif()

if(NOT TARGET dplasma_scalapack_common)
  add_library(dplasma_scalapack_common OBJECT common.c)
  target_include_directories(dplasma_scalapack_common
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(dplasma_scalapack_common
    PUBLIC
    dplasma
    PaRSEC::parsec
    BLACS::blacs)
  set_target_properties(dplasma_scalapack_common PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
  if(SCALAPACK_SUP_UNDERSCORE)
    target_compile_definitions(dplasma_scalapack_common PUBLIC SCALAPACK_SUP_UNDERSCORE)
  endif(SCALAPACK_SUP_UNDERSCORE)
endif(NOT TARGET dplasma_scalapack_common)


### Generate the dplasma wrappers for the scalapack routines
set(SOURCES
  dplasma_wrapper_parsec_init.c
  dplasma_wrapper_pdgemm.c
  dplasma_wrapper_pdgetrf_1d.c
  dplasma_wrapper_pdpotrf.c
  dplasma_wrapper_pdtrmm.c
  dplasma_wrapper_pdtrsm.c
  dplasma_wrapper_pdlatsqr.c
)
# NOT COMPLIANT:
#  dplasma_wrapper_pdgetrf_nopiv.c --> only diagonal dominant matrices

### Generate the lib
if(NOT TARGET dplasma_scalapack)
  if(BUILD_SHARED_LIBS)
    add_library(dplasma_scalapack SHARED ${SOURCES})
  else(BUILD_SHARED_LIBS)
    add_library(dplasma_scalapack STATIC ${SOURCES})
  endif(BUILD_SHARED_LIBS)
endif(NOT TARGET dplasma_scalapack)

set_target_properties(dplasma_scalapack PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

target_link_libraries(dplasma_scalapack
  PRIVATE dplasma_scalapack_common
  PUBLIC
  dplasma
  PaRSEC::parsec)
set_target_properties(dplasma_scalapack PROPERTIES VERSION ${DPLASMA_VERSION_MAJOR}.${DPLASMA_VERSION_MINOR}
                                       SOVERSION ${DPLASMA_VERSION_MAJOR})

install(TARGETS dplasma_scalapack
    EXPORT dplasma_scalapack-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/)
