include(PrecisionGenerator)

set(DTD_HEADERS
  dplasma_z_dtd.h
)
set(DTD_SOURCES
    zgemm.c
    ztrsm.c
    zherk.c
    zpotrf.c
)

set(generated_headers "")
precisions_rules_py(generated_dtd_headers
  "${DTD_HEADERS}"
  PRECISIONS "${DPLASMA_PRECISIONS}")

set(generated_files "")
precisions_rules_py(generated_dtd_cores
  "${DTD_SOURCES}"
  PRECISIONS "${DPLASMA_PRECISIONS}")

add_custom_target(dplasma_dtd_includes ALL SOURCES
                  ${generated_dtd_headers} )

### Publish the documented files
#add_documented_files(PROJECT DPLASMA DIR ${CMAKE_CURRENT_BINARY_DIR} FILES ${generated_dtd_files} ${generated_dtd_headers})

### Generate the dplasma_dtd object library
# We cannot do a simple target_sources because of a bug in CMake <3.18 where the
# GENERATED property has a directory visibility
if(NOT TARGET dplasma_dtd)
  add_library(dplasma_dtd OBJECT ${generated_dtd_cores})
endif(NOT TARGET dplasma_dtd)

add_dependencies(dplasma_dtd dplasma_includes dplasma_cores_includes dplasma_dtd_includes)

set_target_properties(dplasma_dtd PROPERTIES ENABLE_EXPORTS True)
set_target_properties(dplasma_dtd PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
target_include_directories(dplasma_dtd
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
  $<$<NOT:${DPLASMA_BUILD_INPLACE}>:${CMAKE_CURRENT_SOURCE_DIR}>)
  target_link_libraries(dplasma_dtd
    PUBLIC
    dplasma_cores
    PaRSEC::parsec
    LAPACKE::LAPACKE
    $<$<BOOL:${DPLASMA_HAVE_CUDA}>:CUDA::cublas>
    )

# Integrate the output into the main library
target_sources(dplasma PRIVATE $<TARGET_OBJECTS:dplasma_dtd>)