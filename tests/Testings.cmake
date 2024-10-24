#
# Copyright (c) 2010-2024 The University of Tennessee and The University
#                         of Tennessee Research Foundation.  All rights
#                         reserved.
#
#
# A more compact representation of the DPLASMA tests. We can compose any number
# of tests, by providing the epilogue in ALL_TESTS and have the corresponding
# OPTION_** defined to the extra options necessary to run that particular flavor
# of the test.
#
set(OPTIONS_PTG_to_DTD "--;--mca;mca_pins;ptg_to_dtd")
set(DEFAULT_OPTIONS "-x;-v=5")
#set(OPTIONS "")

set(PTG2DTD "ptg2dtd")
set(PTG2DTD_OPTIONS "--;--mca;mca_pins;ptg_to_dtd")
set(OPTIONS "-x;-v=5")
#set(OPTIONS "")

macro(dplasma_add_test m_nameradix m_dependsradix m_types)
  foreach (m_type "${m_types}")
    set(labels "dplasma")
    if (m_type MATCHES "gpu")
      set(m_gpus "${CTEST_GPU_LAUNCHER_OPTIONS}")
      list(APPEND labels "gpu")
    else()
      unset(m_gpus)
    endif()
    if (m_type MATCHES "shm")
      list(APPEND labels "shm")
      string(REPLACE "${PTG2DTD}_" "" m_suffix ${m_type})
      add_test(dplasma_${prec}${m_nameradix}_${m_suffix} ${SHM_TEST_CMD_LIST} ${m_gpus} ./testing_${prec}${m_nameradix} ${ARGN})
      set_tests_properties(dplasma_${prec}${m_nameradix}_${m_suffix} PROPERTIES LABELS "${labels}")
      if (NOT "" STREQUAL "${m_dependsradix}")
        set_tests_properties(dplasma_${prec}${m_nameradix}_${m_suffix} PROPERTIES DEPENDS "launcher_shm;dplasma_${prec}${m_dependsradix}_shm")
      endif()
      if (m_type MATCHES ${PTG2DTD})
        add_test(dplasma_${prec}${m_nameradix}_${PTG2DTD}_${m_suffix} ${SHM_TEST_CMD_LIST} ${m_gpus} ./testing_${prec}${m_nameradix} ${ARGN} ${PTG2DTD_OPTIONS})
        set_tests_properties(dplasma_${prec}${m_nameradix}_${PTG2DTD}_${m_suffix} PROPERTIES DEPENDS dplasma_${prec}${m_nameradix}_${m_suffix} LABELS "${labels};${PTG2DTD}")
      endif()
    endif()

    if (m_type MATCHES "mpi")
      list(APPEND labels "mpi")
      string(REGEX REPLACE ".*mpi:" "" m_procs ${m_type})
      string(REGEX REPLACE ":.*" "" m_suffix ${m_type})
      add_test(dplasma_${prec}${m_nameradix}_${m_suffix} ${MPI_TEST_CMD_LIST} ${m_procs} ${m_gpus} ./testing_${prec}${m_nameradix} ${ARGN})
      set_tests_properties(dplasma_${prec}${m_nameradix}_${m_suffix} PROPERTIES LABELS "${labels}")
      if (NOT "" STREQUAL "${m_dependsradix}")
        set_tests_properties(dplasma_${prec}${m_nameradix}_${m_suffix} PROPERTIES DEPENDS "launcher_mpi;dplasma_${prec}${m_dependsradix}_mpi")
      else()
        set_tests_properties(dplasma_${prec}${m_nameradix}_mpi PROPERTIES DEPENDS launcher_mpi)
      endif()
    endif()
    # enable devices only in tests that explicitely require them
    # restrict memory use for oversubscribed runners
    set_tests_properties(dplasma_${prec}${m_nameradix}_${m_suffix} PROPERTIES ENVIRONMENT
      "PARSEC_MCA_device_cuda_enabled=0;PARSEC_MCA_device_hip_enabled=0;PARSEC_MCA_device_level_zero_enabled=0;PARSEC_MCA_device_cuda_memory_use=10;PARSEC_MCA_device_hip_memory_use=10;PARSEC_MCA_device_level_zero_memory_use=10")
  endforeach()
endmacro()


# The space in the ALL_TESTS list is there to provide room for an empty
# element in the list.
set(ALL_TESTS " ;")
#set(ALL_TESTS " ;_PTG_to_DTD")
#
# Check BLAS/Lapack subroutines in shared memory
#
foreach(prec ${DPLASMA_PRECISIONS} )

  foreach(test ${ALL_TESTS})
    set(OPTIONS "${DEFAULT_OPTIONS};${OPTIONS${test}}")

    # check the control and test matrices generation (zplrnt, zplghe, zplgsy, zpltmg) in shared memory
    dplasma_add_test(print              ""      shm -N 64 -t 7 ${OPTIONS})

    # check the norms that are used in all other testings
    dplasma_add_test(lange              print   shm -M 287 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(lanm2              print   shm -M 287 -N 283 -K 97 -t 56 ${OPTIONS})
    # Need to add testings on zlacpy, zlaset, zgeadd, zlascal, zger, (zlaswp?)

    # BLAS Shared memory
    dplasma_add_test(gemm               lange   shm -M 106 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(symm               lange   shm -M 106 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(trmm               lange   shm -M 106 -N 150 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(trsm               trmm    shm -M 106 -N 150 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(syrk               lange   shm -M 287 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(syr2k              lange   shm -M 287 -N 283 -K 97 -t 56 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(hemm             lange   shm -M 106 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(herk             lange   shm -M 287 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(her2k            lange   shm -M 287 -N 283 -K 97 -t 56 ${OPTIONS})
    endif()

    # Cholesky
    dplasma_add_test(potrf              trsm    shm -N 378 -t 93        ${OPTIONS})
    dplasma_add_test(posv               trsm    shm -N 378 -t 93 -K 367 ${OPTIONS})

    # QR / LQ
    dplasma_add_test(geqrf              gemm    shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(gelqf              gemm    shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(unmqr            gemm    shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(unmlq            gemm    shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    else()
      dplasma_add_test(ormqr            gemm    shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(ormlq            gemm    shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    endif()

    # QR / LQ: HQR
    dplasma_add_test(geqrf_hqr          gemm    shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(gelqf_hqr          gemm    shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(unmqr_hqr        gemm   shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(unmlq_hqr        gemm   shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    else()
      dplasma_add_test(ormqr_hqr        gemm   shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(ormlq_hqr        gemm   shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    endif()

    # QR / LQ: systolic
    dplasma_add_test(geqrf_systolic     gemm shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
    dplasma_add_test(gelqf_systolic     gemm shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(unmqr_systolic   gemm shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(unmlq_systolic   gemm shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    else()
      dplasma_add_test(ormqr_systolic   gemm shm -M 487 -N 283 -K 97 -t 56 ${OPTIONS})
      dplasma_add_test(ormlq_systolic   gemm shm -M 287 -N 383 -K 97 -t 56 ${OPTIONS})
    endif()

    # LU
    dplasma_add_test(getrf_nopiv        gemm    shm -N 378 -t 93       ${OPTIONS})
    dplasma_add_test(getrf_1d           gemm    shm -N 378 -t 93       ${OPTIONS})
    dplasma_add_test(getrf_ptgpanel     gemm    shm -N 378 -t 93       ${OPTIONS})
    dplasma_add_test(getrf_incpiv       gemm    shm -N 378 -t 93 -i 17 ${OPTIONS})
    dplasma_add_test(getrf_qrf          gemm    shm -N 378 -t 93 -i 17 ${OPTIONS})

    #dplasma_add_test(gesv_incpiv        gemm    shm -N 874 -K 367 -t 76       ${OPTIONS})
    dplasma_add_test(gesv_incpiv        gemm    shm -N 874 -K 367 -t 76 -i 23 ${OPTIONS})
  endforeach(test ${ALL_TESTS})

  # Reset the OPTIONS to default values
  set(OPTIONS "${DEFAULT_OPTIONS}")

  # The insert_task interface
  dplasma_add_test(potrf_dtd            trsm    shm -N 874 -K 367 -t 76 -i 23 ${OPTIONS})
  dplasma_add_test(geqrf_dtd            gemm    shm -N 874 -K 367 -t 76 -i 23 ${OPTIONS})
  dplasma_add_test(getrf_incpiv_dtd     gemm    shm -N 874 -K 367 -t 76 -i 23 ${OPTIONS})

  # GPU tests
  if (DPLASMA_HAVE_CUDA)
    dplasma_add_test(potrf              potrf   1gpu_cuda_shm -N 3200 -t 320 ${OPTIONS} -g 1 -- --mca device_cuda_memory_number_of_blocks 4096)
    dplasma_add_test(potrf              potrf   1gpu_cuda_lowmem_shm -N 3200 -t 320 ${OPTIONS} -g 1 -- --mca device_cuda_memory_number_of_blocks 21)
    dplasma_add_test(potrf              potrf   1gpu_cuda_~knb_shm -N 1700 -t 320 ${OPTIONS} -g 1 -- --mca device_cuda_memory_number_of_blocks 4096)
    dplasma_add_test(potrf              potrf   2gpu_cuda_shm -N 4600 -t 320 ${OPTIONS} -g 2 -- --mca device_cuda_memory_number_of_blocks 4096)
    dplasma_add_test(gemm               gemm    1gpu_cuda_shm -N 1280 -t 320 ${OPTIONS} -g 1 -- --mca device_cuda_memory_number_of_blocks 4096)
    dplasma_add_test(gemm               gemm    1gpu_cuda_~knb_shm -N 1000 -t 320 ${OPTIONS} -g 1 -- --mca device_cuda_memory_number_of_blocks 4096)
    dplasma_add_test(gemm               gemm    2gpu_cuda_shm -N 1940 -t 320 ${OPTIONS} -g 2 -- --mca device_cuda_memory_number_of_blocks 4096)
  endif (DPLASMA_HAVE_CUDA)
  if (DPLASMA_HAVE_HIP)
    dplasma_add_test(potrf              potrf   1gpu_hip_shm -N 3200 -t 320 ${OPTIONS} -g 1 -- --mca device_hip_memory_number_of_blocks 4096)
    dplasma_add_test(potrf              potrf   1gpu_hip_lowmem_shm -N 3200 -t 320 ${OPTIONS} -g 1 -- --mca device_hip_memory_number_of_blocks 21)
    dplasma_add_test(potrf              potrf   1gpu_hip_~knb_shm -N 1700 -t 320 ${OPTIONS} -g 1 -- --mca device_hip_memory_number_of_blocks 4096)
    dplasma_add_test(potrf              potrf   2gpu_hip_shm -N 4600 -t 320 ${OPTIONS} -g 2 -- --mca device_hip_memory_number_of_blocks 4096)
    dplasma_add_test(gemm               gemm    1gpu_hip_shm -N 1280 -t 320 ${OPTIONS} -g 1 -- --mca device_hip_memory_number_of_blocks 4096)
    dplasma_add_test(gemm               gemm    1gpu_hip_~knb_shm -N 1000 -t 320 ${OPTIONS} -g 1 -- --mca device_hip_memory_number_of_blocks 4096)
    dplasma_add_test(gemm               gemm    2gpu_hip_shm -N 1940 -t 320 ${OPTIONS} -g 2 -- --mca device_hip_memory_number_of_blocks 4096)
  endif (DPLASMA_HAVE_HIP)


  #   if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
  #     dplasma_add_test(heev "" ${PTG2DTD}_shm -N 4000 ${OPTIONS})
  #   else()
  #     dplasma_add_test(syev "" ${PTG2DTD}_shm -N 4000 ${OPTIONS})
  #   endif()
  #   dplasma_add_test(geqrf_p0 "" ${PTG2DTD}_shm -N 4000 -t 200 -i 32 -x --qr_a=2 --treel 0 --tsrr=0 -v=5)
  #   dplasma_add_test(geqrf_p1 "" ${PTG2DTD}_shm -N 4000 -t 200 -i 32 -x --qr_a=2 --treel 1 --tsrr=0 -v=5)
  #   dplasma_add_test(geqrf_p2 "" ${PTG2DTD}_shm -N 4000 -t 200 -i 32 -x --qr_a=2 --treel 2 --tsrr=0 -v=5)
  #   dplasma_add_test(geqrf_p3 "" ${PTG2DTD}_shm -N 4000 -t 200 -i 32 -x --qr_a=2 --treel 3 --tsrr=0 -v=5)
  #
endforeach()

#
# Distributed Memory Testings
#
if( MPI_C_FOUND )
  set(PROCS 4)
  set(CORES "")
  #set(CORES "-c;1")

  foreach(prec ${DPLASMA_PRECISIONS})

    # check the control and test matrices generation (zplrnt, zplghe, zplgsy, zpltmg) in distributed memory
    dplasma_add_test(print              ""    mpi:${PROCS} -N 64 -t 7 ${OPTIONS})

    # check the norms that are used in all other testings
    dplasma_add_test(lange              print mpi:${PROCS} -M 287 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(lanm2              print mpi:${PROCS} -M 287 -N 283 -K 97 -t 19 ${OPTIONS})

    # Need to add testings on zlacpy, zlaset, zgeadd, zlascal, zger, (zlaswp?)

    # BLAS Shared memory
    dplasma_add_test(trmm               lange mpi:${PROCS} -M 106 -N 150 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(trsm               trmm  mpi:${PROCS} -M 106 -N 150 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(gemm               lange mpi:${PROCS} -M 106 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(gemm_dtd           lange mpi:${PROCS} -M 106 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(symm               lange mpi:${PROCS} -M 106 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(syrk               lange mpi:${PROCS} -M 287 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(syr2k              lange mpi:${PROCS} -M 287 -N 283 -K 97 -t 19 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(hemm             lange mpi:${PROCS} -M 106 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(herk             lange mpi:${PROCS} -M 287 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(her2k            lange mpi:${PROCS} -M 287 -N 283 -K 97 -t 19 ${OPTIONS})
    endif()

    # Cholesky
    dplasma_add_test(potrf              gemm mpi:${PROCS} -N 378 -t 19        ${OPTIONS})
    dplasma_add_test(potrf_dtd          gemm mpi:${PROCS} -N 378 -t 19        ${OPTIONS})
    dplasma_add_test(potrf_dtd_untied   gemm mpi:${PROCS} -N 378 -t 19        ${OPTIONS})
    dplasma_add_test(posv               gemm mpi:${PROCS} -N 457 -t 19 -K 367 ${OPTIONS})

    # QR / LQ
    dplasma_add_test(geqrf              gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(geqrf_dtd          gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(gelqf              gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(unmqr            gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(unmlq            gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    else()
      dplasma_add_test(ormqr            gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(ormlq            gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    endif()

    # QR / LQ: HQR
    dplasma_add_test(geqrf_hqr          gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(gelqf_hqr          gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(unmqr_hqr        gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(unmlq_hqr        gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    else()
      dplasma_add_test(ormqr_hqr        gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(ormlq_hqr        gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    endif()

    # QR / LQ: systolic
    dplasma_add_test(geqrf_systolic     gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
    dplasma_add_test(gelqf_systolic     gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
      dplasma_add_test(unmqr_systolic   gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(unmlq_systolic   gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    else()
      dplasma_add_test(ormqr_systolic   gemm mpi:${PROCS} -M 487 -N 283 -K 97 -t 19 ${OPTIONS})
      dplasma_add_test(ormlq_systolic   gemm mpi:${PROCS} -M 287 -N 383 -K 97 -t 19 ${OPTIONS})
    endif()

    # LU
    dplasma_add_test(getrf_1d           gemm mpi:${PROCS} -N 378 -t 19 -P 1 ${OPTIONS})
    dplasma_add_test(getrf_incpiv       gemm mpi:${PROCS} -N 378 -t 19 -i 7 ${OPTIONS})
    dplasma_add_test(getrf_incpiv_dtd   gemm mpi:${PROCS} -N 378 -t 19 -i 7 ${OPTIONS})
    dplasma_add_test(getrf_ptgpanel     gemm mpi:${PROCS} -N 378 -t 19      ${OPTIONS})
    dplasma_add_test(getrf_nopiv        gemm mpi:${PROCS} -N 378 -t 19      ${OPTIONS})
    dplasma_add_test(getrf_qrf          gemm mpi:${PROCS} -N 378 -t 19 -i 7 ${OPTIONS})

    #dplasma_add_test(gesv "" mpi:${PROCS} -N 874 -K 367 -t 76       ${OPTIONS})
    dplasma_add_test(gesv_incpiv        gemm mpi:${PROCS} -N 874 -K 367 -t 17 -i 7 ${OPTIONS})

    # GPU Cholesky tests
    if (DPLASMA_HAVE_CUDA AND MPI_C_FOUND)
        dplasma_add_test(potrf potrf      1gpu_cuda_mpi:${PROCS} -N 3200 -t 320 ${OPTIONS} -g 1 -P 2 -- --mca device_cuda_memory_number_of_blocks 4096)
        dplasma_add_test(potrf potrf      1gpu_cuda_~knb_mpi:${PROCS} -N 1700 -t 320 ${OPTIONS} -g 1 -P 2 -- --mca device_cuda_memory_number_of_blocks 4096)
        dplasma_add_test(potrf potrf_1gpu 2gpu_cuda_mpi:${PROCS} -N 4600 -t 320 ${OPTIONS} -g 2 -P 2 -- --mca device_cuda_memory_number_of_blocks 4096)
        dplasma_add_test(gemm  gemm       2gpu_cuda_mpi:${PROCS} -N 1940 -t 320 ${OPTIONS} -g 2 -P 2 -- --mca device_cuda_memory_number_of_blocks 4096)
        dplasma_add_test(gemm  gemm       2gpu_cuda_lowmem_mpi:${PROCS} -N 1940 -t 320 ${OPTIONS} -g 2 -P 2 -- --mca device_cuda_memory_number_of_blocks 21)
    endif (DPLASMA_HAVE_CUDA AND MPI_C_FOUND)
    if (DPLASMA_HAVE_HIP AND MPI_C_FOUND)
        dplasma_add_test(potrf potrf      1gpu_hip_mpi:${PROCS} -N 3200 -t 320 ${OPTIONS} -g 1 -P 2 -- --mca device_hip_memory_number_of_blocks 4096)
        dplasma_add_test(potrf potrf      1gpu_hip_~knb_mpi:${PROCS} -N 1700 -t 320 ${OPTIONS} -g 1 -P 2 -- --mca device_hip_memory_number_of_blocks 4096)
        dplasma_add_test(potrf potrf_1gpu 2gpu_hip_mpi:${PROCS} -N 4600 -t 320 ${OPTIONS} -g 2 -P 2 -- --mca device_hip_memory_number_of_blocks 4096)
        dplasma_add_test(gemm  gemm       2gpu_hip_mpi:${PROCS} -N 1940 -t 320 ${OPTIONS} -g 2 -P 2 -- --mca device_hip_memory_number_of_blocks 4096)
        dplasma_add_test(gemm  gemm       2gpu_hip_lowmem_mpi:${PROCS} -N 1940 -t 320 ${OPTIONS} -g 2 -P 2 -- --mca device_hip_memory_number_of_blocks 21)
    endif (DPLASMA_HAVE_HIP AND MPI_C_FOUND)

    # dplasma_add_test(potrf_pbq "" mpi:${PROCS} -N 4000 ${OPTIONS} -o PBQ)
    # dplasma_add_test(geqrf_pbq "" mpi:${PROCS} -N 4000 ${OPTIONS} -o PBQ)
    # dplasma_add_test(geqrf_p0 "" mpi:${PROCS} -N 4000 -t 200 -i 32 -x --qr_p=4 --qr_a=2 --treel 0 --tsrr=0 -v=5)
    # dplasma_add_test(geqrf_p1 "" mpi:${PROCS} -N 4000 -t 200 -i 32 -x --qr_p=4 --qr_a=2 --treel 1 --tsrr=0 -v=5)
    # dplasma_add_test(geqrf_p2 "" mpi:${PROCS} -N 4000 -t 200 -i 32 -x --qr_p=4 --qr_a=2 --treel 2 --tsrr=0 -v=5)
    # dplasma_add_test(geqrf_p3 "" mpi:${PROCS} -N 4000 -t 200 -i 32 -x --qr_p=4 --qr_a=2 --treel 3 --tsrr=0 -v=5)

    # if ( ${prec} STREQUAL "c" OR ${prec} STREQUAL "z" )
    #   dplasma_add_test(heev "" mpi:${PROCS} -N 2000 ${OPTIONS})
    # else()
    #   dplasma_add_test(syev "" mpi:${PROCS} -N 2000 ${OPTIONS})
    # endif()
  endforeach()

endif( MPI_C_FOUND )

