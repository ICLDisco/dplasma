include(PrecisionGenerator)
# reset variables
set(generated_files "")
set(generated_headers "")
set(generated_wrappers "")
set(cplx_files "")
set(all_precisions_files "")

set(CORE_CUDA_HEADERS "")
if( CUDA_FOUND )
    set(CORE_CUDA_SOURCES dplasma_cuda_ztsmqr.c)
else()
    set(CORE_CUDA_SOURCES "")
endif()

### generate the dplasma_cores headers for all required precisions
set(HEADERS
    dplasma_zcores.h
    ${CORE_CUDA_HEADERS})
precisions_rules_py(generated_headers
                    ${HEADERS}
                    PRECISIONS "s;d;c;z")
add_custom_target(dplasma_cores_includes ALL SOURCES
                  ${generated_headers} )

foreach(generated_header ${generated_headers})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${generated_header}
          DESTINATION ${DPLASMA_INSTALL_PREFIX}/include/cores)
endforeach()

### Generate the dplasma wrappers for all required precisions
set(SOURCES
  core_dzasum.c
  core_zamax.c
  core_zgeadd.c
  core_zgelqt.c
  core_zgemdm.c
  core_zgemm.c
  core_zgemv.c
  core_zgeqrt.c
  core_zgessm.c
  core_zgessq.c
  core_zgetf2_nopiv.c
  core_zgetrf_incpiv.c
  core_zgetrf_nopiv.c
  core_zgetrf_reclap.c
  core_zgetrf_rectil.c
  core_zhebut.c
  core_zhedrk.c
  core_zhemm.c
  core_zher2k.c
  core_zherfb.c
  core_zherk.c
  core_zhessq.c
  core_zhetrf2_nopiv.c
  core_zhetrf_nopiv.c
  core_zlacpy.c
  core_zlange.c
  core_zlanhe.c
  core_zlansy.c
  core_zlantr.c
  core_zlaset.c
  core_zlaswp.c
  core_zlauum.c
  core_zpamm.c
  core_zparfb.c
  core_zpemv.c
  core_zplghe.c
  core_zplgsy.c
  core_zplrnt.c
  core_zpltmg.c
  core_zpltmg_chebvand.c
  core_zpltmg_circul.c
  core_zpltmg_condex.c
  core_zpltmg_fiedler.c
  core_zpltmg_hankel.c
  core_zpltmg_toeppd.c
  core_zpotrf.c
  core_zshift.c
  core_zssssm.c
  core_zsymm.c
  core_zsyr2k.c
  core_zsyrk.c
  core_zsyssq.c
  core_ztradd.c
  core_ztrasm.c
  core_ztrdv.c
  core_ztrmdm.c
  core_ztrmm.c
  core_ztrsm.c
  core_ztrssq.c
  core_ztrtri.c
  core_ztslqt.c
  core_ztsmlq.c
  core_ztsmlq_hetra1.c
  core_ztsmqr.c
  core_ztsmqr_corner.c
  core_ztsmqr_hetra1.c
  core_ztsqrt.c
  core_ztstrf.c
  core_zttlqt.c
  core_zttmlq.c
  core_zttmqr.c
  core_zttqrt.c
  core_zunmlq.c
  core_zunmqr.c
)
set(ALL_PREC_SRCS core_zgessq.c )

set(CPLX_ONLY core_dlag2z.c )

precisions_rules_py(generated_files
    "${SOURCES}" PRECISIONS "${DPLASMA_PRECISIONS}")

precisions_rules_py(all_precisions_files
    "${ALL_PREC_SRCS}" PRECISIONS "s;d;c;z")

precisions_rules_py(cplx_files
    "${CPLX_ONLY}" PRECISIONS "c;z")

set(common_srcs
    gkkleader.c
    primes.c
    descriptor.c
    global.c)

### Generate the dplasma_cores archive
if(NOT TARGET dplasma_cores_obj)
  add_library(dplasma_cores_obj OBJECT ${generated_files} ${all_precisions_files} ${cplx_files} ${common_srcs})
endif(NOT TARGET dplasma_cores_obj)

target_include_directories(dplasma_cores_obj
  PRIVATE
  ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
  $<$<NOT:${DPLASMA_BUILD_INPLACE}>:${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/include>)

set_target_properties(dplasma_cores_obj PROPERTIES ENABLE_EXPORTS True)
set_target_properties(dplasma_cores_obj PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

add_dependencies(dplasma_cores_obj dplasma_includes dplasma_cores_includes)

add_documented_files(PROJECT DPLASMA DIR ${CMAKE_CURRENT_BINARY_DIR} FILES ${generated_files} ${all_precisions_files} ${cplx_files})
add_documented_files(PROJECT DPLASMA DIR ${CMAKE_CURRENT_SOURCE_DIR} FILES ${common_srcs})

if( CUDA_FOUND )
  # generate the cores_cuda library (hooks to cublas)
  precisions_rules_py(generated_cuda_files
                      ${CORE_CUDA_SOURCES}
                      PRECISIONS "${DPLASMA_PRECISIONS}")

  add_library(dplasma_cores_cuda OBJECT ${generated_cuda_files})
  add_dependencies(dplasma_cores_cuda dplasma_includes
                                      dplasma_cores_includes)
  target_include_directories(dplasma_cores_cuda
    PUBLIC
    ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
    $<$<NOT:${DPLASMA_BUILD_INPLACE}>:${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/include>)

  ### Generate the CUDA .cu kernels if necessary
  # By default drop support for older devices (add 11 and 13 for support
  # SM targets with poor DP flops are also disabled by default (add 30 50 52 53 61 62)
  if( ${CUDA_VERSION_STRING} VERSION_LESS "4.2" )
    set(CUDA_SM_TARGETS 11 13 20 CACHE STRING "Specify GPU architectures to build binaries for, BIN(PTX) format is supported")
elseif( ${CUDA_VERSION_STRING} VERSION_LESS "8.0" )
    set(CUDA_SM_TARGETS 20 35 37 CACHE STRING "Specify GPU architectures to build binaries for, BIN(PTX) format is supported")
  elseif( ${CUDA_VERSION_STRING} VERSION_EQUAL "7.0" )
    set(CUDA_SM_TARGETS 20 35 37 60 CACHE STRING "Specify GPU architectures to build binaries for, BIN(PTX) format is supported")
  else()
    # SM_20 will be removed in cuda>8
    set(CUDA_SM_TARGETS 35 37 60 CACHE STRING "Specify GPU architectures to build binaries for, BIN(PTX) format is supported")
  endif()

# For Tesla
  set(CUDA_NVCC_FLAGS_11 -maxrregcount 32 -arch sm_11 -DCUDA_SM_VERSION=11)
  set(CUDA_NVCC_FLAGS_13 -maxrregcount 32 -arch sm_13 -DCUDA_SM_VERSION=13)
# For Fermi
  set(CUDA_NVCC_FLAGS_20                  -arch sm_20 -DCUDA_SM_VERSION=20)
# For Kepler
  set(CUDA_NVCC_FLAGS_30                  -arch sm_30 -DCUDA_SM_VERSION=30)
  set(CUDA_NVCC_FLAGS_35                  -arch sm_35 -DCUDA_SM_VERSION=35)
  set(CUDA_NVCC_FLAGS_37                  -arch sm_37 -DCUDA_SM_VERSION=37)
# For Maxwell
  set(CUDA_NVCC_FLAGS_50                  -arch sm_50 -DCUDA_SM_VERSION=50)
  set(CUDA_NVCC_FLAGS_52                  -arch sm_52 -DCUDA_SM_VERSION=52)
  set(CUDA_NVCC_FLAGS_53                  -arch sm_53 -DCUDA_SM_VERSION=53)
# For Pascal
  set(CUDA_NVCC_FLAGS_60                  -arch sm_60 -DCUDA_SM_VERSION=60)
  set(CUDA_NVCC_FLAGS_61                  -arch sm_61 -DCUDA_SM_VERSION=61)
  set(CUDA_NVCC_FLAGS_62                  -arch sm_62 -DCUDA_SM_VERSION=62)

  ### CUDA .cu sources
  set( CUDA11_SOURCES
#       zgemm_11_12_13.cu
  )
  set( CUDA13_SOURCES ${CUDA11_SOURCES} )
  set( CUDA20_SOURCES
#        zgemm_20_30.cu
  )
  set( CUDA30_SOURCES ${CUDA20_SOURCES} )
  set( CUDA35_SOURCES ${CUDA30_SOURCES} )
  set( CUDA37_SOURCES ${CUDA30_SOURCES} )
  set( CUDA50_SOURCES ${CUDA30_SOURCES} )
  set( CUDA52_SOURCES ${CUDA50_SOURCES} )
  set( CUDA53_SOURCES ${CUDA50_SOURCES} )
  set( CUDA60_SOURCES ${CUDA50_SOURCES} )
  set( CUDA61_SOURCES ${CUDA60_SOURCES} )
  set( CUDA62_SOURCES ${CUDA60_SOURCES} )

  foreach( _smtarget ${CUDA_SM_TARGETS} )
    if( NOT "${CUDA${_smtarget}_SOURCES}" STREQUAL "" )
      set(cuda_generated_files "")
      precisions_rules_py(cuda_generated_files
                          "${CUDA${_smtarget}_SOURCES}"
                          PRECISIONS "${DPLASMA_PRECISIONS}")
      foreach( _cudafile ${cuda_generated_files} )
          get_source_file_property(_IsInBinaryDir ${_cudafile} IS_IN_BINARY_DIR )
          if( _IsInBinaryDir )
              string(REGEX REPLACE
                      ${_cudafile}
                      ${CMAKE_CURRENT_BINARY_DIR}/${_cudafile}
                      cuda_generated_files "${cuda_generated_files}")
          endif( _IsInBinaryDir )
      endforeach()

      set(CUDA_NVCC_FLAGS_BACKUP ${CUDA_NVCC_FLAGS})
      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS_${_smtarget}} ${CUDA_NVCC_FLAGS})
      cuda_add_library(dplasma_cucores_sm${_smtarget}
                       ${cuda_generated_files})
      cuda_add_cublas_to_target( dplasma_cucores_sm${_smtarget} )
      install(TARGETS dplasma_cucores_sm${_smtarget}
        LIBRARY DESTINATION ${DPLASMA_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${DPLASMA_INSTALL_PREFIX}/lib )
      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS_BACKUP})
      if( NOT CMAKE_BUILD_SHARED_LIBS )
          target_link_libraries( dplasma_cores_cuda dplasma_cucores_sm${_smtarget} )
      endif()
    endif()
  endforeach()

  set_target_properties(dplasma_cores_cuda PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
  cuda_build_clean_target()
endif()

if(NOT TARGET dplasa_cores)
  add_library(dplasma_cores INTERFACE)
endif(NOT TARGET dplasa_cores)

target_include_directories(dplasma_cores
  INTERFACE
  ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
  $<$<NOT:${DPLASMA_BUILD_INPLACE}>:${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/include>)

target_sources(dplasma_cores INTERFACE
               $<TARGET_OBJECTS:dplasma_cores_obj>
               $<$<BOOL:${DPLASMA_HAVE_CUDA}>:$<TARGET_OBJECTS:dplasma_cores_cuda>>)


target_link_libraries(dplasma_cores INTERFACE
  m
  PaRSEC::parsec
  ${LAPACKE_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES})
if( ${DPLASMA_HAVE_CUDA} )
  target_link_libraries(dplasma_cores
      INTERFACE ${CUDA_CUBLAS_LIBRARIES})
endif( ${DPLASMA_HAVE_CUDA} )
add_documented_files(PROJECT DPLASMA DIR ${CMAKE_CURRENT_BINARY_DIR} FILES ${generated_files} ${all_precisions_files} ${cplx_files})
add_documented_files(PROJECT DPLASMA DIR ${CMAKE_CURRENT_SOURCE_DIR} FILES ${common_srcs})

