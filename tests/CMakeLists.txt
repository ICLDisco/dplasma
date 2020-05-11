include(DplasmaCompilerFlags)
if(NOT TARGET dplasma_tests_common)
  add_library(dplasma_tests_common OBJECT common.c)
  target_include_directories(dplasma_tests_common
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
#  common_timing.h flops.h)
  target_link_libraries(dplasma_tests_common
    PUBLIC
    dplasma
    $<$<BOOL:${MPI_C_FOUND}>:MPI::MPI_C>)
  add_dependencies(dplasma_tests_common dplasma_includes)
endif(NOT TARGET dplasma_tests_common)

include(RulesTestings.cmake)

set(TESTS
   # Lapack Auxiliary
   testing_zgeadd.c
   testing_zprint.c
   testing_zlange.c
   testing_zlanm2.c
   # Level 3 Blas
   testing_zgemm.c
   testing_zhemm.c
   testing_zsymm.c
   testing_zherk.c
   testing_zsyrk.c
   testing_zher2k.c
   testing_zsyr2k.c
   testing_ztrmm.c
   testing_ztrsm.c
   #
   # Lapack
   #
   # Cholesky
   testing_zpotrf.c
   testing_zpotrf_fullm.c
   testing_zposv.c
   testing_ztrtri.c
   testing_zpoinv.c
   # QR/LQ
   testing_zgelqf.c
   testing_zgeqrf.c
   testing_zunmqr.c
   testing_zunmlq.c
   testing_zgeqrf_hqr.c
   testing_zgelqf_hqr.c
   testing_zunmqr_hqr.c
   testing_zunmlq_hqr.c
   testing_zgeqrf_systolic.c
   testing_zgelqf_systolic.c
   testing_zunmqr_systolic.c
   testing_zunmlq_systolic.c
   # LU
   testing_zgetrf.c
   testing_zgetrf_incpiv.c
   testing_zgetrf_nopiv.c
   testing_zgetrf_ptgpanel.c
   testing_zgetrf_ptgpanel_panel.c
   testing_zgetrf_qrf.c
   #testing_zgesv.c
   testing_zgesv_incpiv.c
   # Extra testings
   testing_zgesvd.c
   testing_zgebrd_ge2gb.c
   testing_zgeqrf_rd.c
   testing_zhebut.c
   testing_zhbrdt.c
   testing_zheev.c
   # insert_task interface
   testing_zpotrf_dtd.c
   testing_zpotrf_dtd_untied.c
   testing_zgeqrf_dtd.c
   testing_zgeqrf_dtd_untied.c
   testing_zgetrf_incpiv_dtd.c
   testing_zgemm_dtd.c
) 
if( DPLASMA_DEBUG_QR_PIVGEN )
  list(APPEND TESTS testing_zpivgen.c)
  include(TestsQRPivgen.cmake)
endif()

testings_addexec(targets "${DPLASMA_PRECISIONS}" "${TESTS}")

option(DPLASMA_INSTALL_TESTS OFF "Install the DPLASMA tests (in libexec)")
if(DPLASMA_INSTALL_TESTS)
  install(TARGETS ${targets}
          DESTINATION ${CMAKE_INSTALL_LIBEXECDIR}/dplasma/tests)
endif(DPLASMA_INSTALL_TESTS)

# Define Testings
include(Testings.cmake)
