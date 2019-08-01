include(PrecisionGenerator)

# reset variables
set(generated_headers "")

### generate the dplasma headers for all required precisions
precisions_rules_py(generated_headers
  "dplasma_z.h"
  PRECISIONS "${DPLASMA_PRECISIONS}")

set(DPLASMA_INCLUDE_STRING "")
foreach(generated_header ${generated_headers})
  set(DPLASMA_INCLUDE_STRING "${DPLASMA_INCLUDE_STRING}\n
  #include \"${generated_header}\"")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${generated_header}
          DESTINATION ${DPLASMA_INSTALL_PREFIX}/include)
endforeach()

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/dplasma.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dplasma.h" @ONLY)

if(NOT TARGET dplasma_includes)
   add_library(dplasma_includes INTERFACE)
endif(NOT TARGET dplasma_includes)
 
target_include_directories(dplasma_includes INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(dplasma_includes_gen ALL
                 DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/dplasma.h
                 SOURCES ${generated_headers} )

add_dependencies(dplasma_includes dplasma_includes_gen)
	       
### install the dplasma headers
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/dplasma.h
  dplasma_complex.h
  dplasma_qr_param.h
  dplasmatypes.h
  DESTINATION ${DPLASMA_INSTALL_PREFIX}/include)

add_documented_files(PROJECT "DPLASMA" DIR ${CMAKE_CURRENT_BINARY_DIR} FILES ${generated_headers} "dplasma.h")
add_documented_files(PROJECT "DPLASMA" FILES "dplasma_complex.h" "dplasma_qr_param.h" "dplasmatypes.h")
