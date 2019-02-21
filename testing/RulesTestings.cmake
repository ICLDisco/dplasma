include(PrecisionGenerator)

macro(testings_addexec OUTPUTLIST PRECISIONS ZSOURCES)
  set(testings_addexec_GENFILES "")
  precisions_rules_py(testings_addexec_GENFILES
    "${ZSOURCES}"
    PRECISIONS "${PRECISIONS}")
  foreach(testings_addexec_GENFILE ${testings_addexec_GENFILES})
    string(REGEX REPLACE "\\.c" "" testings_addexec_EXEC ${testings_addexec_GENFILE})

    if(NOT TARGET ${testings_addexec_EXEC})
      add_executable(${testings_addexec_EXEC} ${testings_addexec_GENFILE})
      add_dependencies(${testings_addexec_EXEC} dplasma dplasma_includes)
      target_link_libraries(${testings_addexec_EXEC}
	PRIVATE common dplasma
	$<$<BOOL:${MPI_C_FOUND}>:MPI::MPI_C>)
    endif(NOT TARGET ${testings_addexec_EXEC})

    if( PLASMA_F_COMPILE_SUCCESS )
      set_target_properties(${testings_addexec_EXEC} PROPERTIES
                              LINKER_LANGUAGE Fortran
                              COMPILE_FLAGS "${testings_addexec_CFLAGS}"
                              LINK_FLAGS "${testings_addexec_LDFLAGS} ${LOCAL_FORTRAN_LINK_FLAGS}")
    else( PLASMA_F_COMPILE_SUCCESS )
      set_target_properties(${testings_addexec_EXEC} PROPERTIES
                              COMPILE_FLAGS "${testings_addexec_CFLAGS}"
                              LINK_FLAGS "${testings_addexec_LDFLAGS}")
    endif( PLASMA_F_COMPILE_SUCCESS )

    list(APPEND ${OUTPUTLIST} ${testings_addexec_EXEC})
  endforeach()

endmacro(testings_addexec)

