include (CheckCCompilerFlag)
include (CheckFortranCompilerFlag)

#
# Check compiler debug flags and capabilities
#

# add gdb symbols in debug and relwithdebinfo, g3 for macro support when available
check_c_compiler_flag( "-g3" DPLASMA_HAVE_G3 )
if( DPLASMA_HAVE_G3 )
  set(wflags "-g3")
else()
  set(wflags "-g")
endif()

# Some compilers produce better debugging outputs with Og vs O0
check_c_compiler_flag( "-Og" DPLASMA_HAVE_Og )
if( DPLASMA_HAVE_Og )
  set(o0flag "-Og")
else()
  set(o0flag "-O0")
endif()

# Set warnings for debug builds
check_c_compiler_flag( "-Wall" DPLASMA_HAVE_WALL )
if( DPLASMA_HAVE_WALL )
  list(APPEND wflags "-Wall" )
endif( DPLASMA_HAVE_WALL )
check_c_compiler_flag( "-Wextra" DPLASMA_HAVE_WEXTRA )
if( DPLASMA_HAVE_WEXTRA )
  list(APPEND wflags "-Wextra" )
endif( DPLASMA_HAVE_WEXTRA )

#
# flags for Intel icc
#
string(REGEX MATCH ".*icc$" _match_icc ${CMAKE_C_COMPILER})
if(_match_icc)
  # Silence annoying warnings
  check_c_compiler_flag( "-wd424" DPLASMA_HAVE_WD )
  if( DPLASMA_HAVE_WD )
    # 424: checks for duplicate ";"
    # 981: every volatile triggers a "unspecified evaluation order", obnoxious
    #      but might be useful for some debugging sessions.
    # 1419: warning about extern functions being declared in .c
    #       files
    # 1572: cuda compares floats with 0.0f.
    # 11074: obnoxious about not inlining long functions.
    list(APPEND wflags "-wd424,981,1419,1572,10237,11074,11076")
  endif( DPLASMA_HAVE_WD )
else(_match_icc)
  check_c_compiler_flag( "-Wno-parentheses-equality" DPLASMA_HAVE_PAR_EQUALITY )
  if( DPLASMA_HAVE_PAR_EQUALITY )
    list(APPEND wflags "-Wno-parentheses-equality")
  endif( DPLASMA_HAVE_PAR_EQUALITY )
endif(_match_icc)

# verbose compilation in debug
add_compile_options(
  "$<$<CONFIG:DEBUG>:${o0flag};${wflags}>"
  "$<$<CONFIG:RELWITHDEBINFO>:${wflags}>")
# remove asserts in release
add_compile_definitions(
  $<$<CONFIG:RELEASE>:NDEBUG>)

#
# Fortran tricks: fortran main avoidance
#
IF (CMAKE_Fortran_COMPILER_WORKS)
  STRING(REGEX MATCH "Intel" _match_intel ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_intel)
    STRING(APPEND CMAKE_Fortran_LINK_EXECUTABLE " -nofor-main")
  ENDIF (_match_intel)

  STRING(REGEX MATCH "PGI$" _match_pgi ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_pgi)
    STRING(APPEND CMAKE_Fortran_LINK_EXECUTABLE " -Mnomain -Bstatic")
  ENDIF (_match_pgi)

  STRING(REGEX MATCH ".*xlf$" _match_xlf ${CMAKE_Fortran_COMPILER})
  IF (_match_xlf)
    MESSAGE(ERROR "Please use the thread-safe version of the xlf compiler (xlf_r)")
  ENDIF (_match_xlf)
ENDIF (CMAKE_Fortran_COMPILER_WORKS)
