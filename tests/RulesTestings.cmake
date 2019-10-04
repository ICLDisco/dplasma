# Copyright (c) 2009-2019 The University of Tennessee and The University
#                         of Tennessee Research Foundation.  All rights
#                         reserved.
#
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#

include(PrecisionGenerator)

# Create a list of targets that correspond to the
# files provided in $ARGN for the precisions (sdcz)
# provided in PRECISIONS.
# The created targets are listed in OUTPUTLIST
function(testings_addexec OUTPUTLIST PRECISIONS)
  set(generated_testings "")
  precisions_rules_py(generated_testings
    ${ARGN}
    PRECISIONS "${PRECISIONS}")
  foreach(generated_testing ${generated_testings})
    string(REGEX REPLACE "\\.c" "" exec ${generated_testing})

    if(NOT TARGET ${exec})
      add_executable(${exec} ${generated_testing})
      add_dependencies(${exec} dplasma dplasma_includes)
      target_link_libraries(${exec}
                            PRIVATE common
                            PUBLIC dplasma
                            $<$<BOOL:${MPI_C_FOUND}>:MPI::MPI_C>)
    endif(NOT TARGET ${exec})

    if( CMAKE_Fortran_COMPILER_WORKS )
      set_target_properties(${exec} PROPERTIES
                            LINKER_LANGUAGE Fortran)
    endif()

    list(APPEND outputlist ${exec})
  endforeach()

  set(${OUTPUTLIST} ${outputlist} PARENT_SCOPE)
endfunction(testings_addexec)

