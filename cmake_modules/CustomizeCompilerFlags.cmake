
include (CMakeDetermineSystem)
include (CheckCCompilerFlag)
include (CheckCSourceCompiles)
include (CheckFortranCompilerFlag)
include (CheckFortranSourceCompiles)
include (CheckFunctionExists)
include (CheckSymbolExists)
include (CheckIncludeFiles)
include (CMakePushCheckState)

#
# Fix the building system for 32 or 64 bits.
#
# On MAC OS X there is a easy solution, by setting the
# CMAKE_OSX_ARCHITECTURES to a subset of the following values:
# ppc;ppc64;i386;x86_64.
# On Linux this is a little bit tricky. We have to check that the
# compiler supports the -m32/-m64 flags as well as the linker.
# Once this issue is resolved the CMAKE_C_FLAGS and CMAKE_EXE_LINKER_FLAGS
# have to be updated accordingly.
#
# TODO: Same trick for the Fortran compiler...
#       no idea how to correctly detect if the required/optional
#          libraries are in the correct format.
#
STRING(REGEX MATCH ".*xlc$" _match_xlc ${CMAKE_C_COMPILER})
IF (_match_xlc)
  MESSAGE(ERROR "Please use the thread-safe version of the xlc compiler (xlc_r)")
ENDIF (_match_xlc)
STRING(REGEX MATCH "XL" _match_xlc ${CMAKE_C_COMPILER_ID})
if (BUILD_64bits)
  if( _match_xlc)
    set( ARCH_BUILD "-q64" )
  else (_match_xlc)
    if( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "sparc64fx" )
      set ( ARCH_BUILD " " )
    else()
      set( ARCH_BUILD "-m64" )
    endif()
  endif(_match_xlc)
else (BUILD_64bits)
  if( _match_xlc)
    set( ARCH_BUILD "-q32" )
  else (_match_xlc)
    set( ARCH_BUILD "-m32" )
  endif(_match_xlc)
endif (BUILD_64bits)

check_c_compiler_flag( ${ARCH_BUILD} C_M32or64 )
if( C_M32or64 )
  string(APPEND CMAKE_C_FLAGS " ${ARCH_BUILD}")
  # Try the same for Fortran and CXX:
  # Use the same 64bit flag as the C compiler if possible
  if(CMAKE_Fortran_COMPILER_WORKS)
    check_fortran_compiler_flag( ${ARCH_BUILD} F_M32or64 )
    if( F_M32or64 )
      string(APPEND CMAKE_Fortran_FLAGS " ${ARCH_BUILD}")
    endif( F_M32or64 )
  endif()
  if(CMAKE_CXX_COMPILER_WORKS)
    check_fortran_compiler_flag( ${ARCH_BUILD} CXX_M32or64 )
    if( CXX_M32or64 )
      string(APPEND CMAKE_CXX_FLAGS " ${ARCH_BUILD}")
    endif( CXX_M32or64 )
  endif()
endif( C_M32or64 )

#
# Check compiler flags and capabilities
#

# Set warnings for debug builds
CHECK_C_COMPILER_FLAG( "-Wall" CC_HAS_WALL )
IF( CC_HAS_WALL )
  STRING( APPEND DPLASMA_C_WFLAGS " -Wall" )
ENDIF( CC_HAS_WALL )
CHECK_C_COMPILER_FLAG( "-Wextra" CC_HAS_WEXTRA )
IF( CC_HAS_WEXTRA )
  STRING( APPEND DPLASMA_C_WFLAGS " -Wextra" )
ENDIF( CC_HAS_WEXTRA )

#
# flags for Intel icc
#
STRING(REGEX MATCH ".*icc$" _match_icc ${CMAKE_C_COMPILER})
if (_match_icc)
  # Silence annoying warnings
  CHECK_C_COMPILER_FLAG( "-wd424" CC_HAS_WD )
  IF( CC_HAS_WD )
    # 424: checks for duplicate ";"
    # 981: every volatile triggers a "unspecified evaluation order", obnoxious
    #      but might be useful for some debugging sessions.
    # 1419: warning about extern functions being declared in .c
    #       files
    # 1572: cuda compares floats with 0.0f.
    # 11074: obnoxious about not inlining long functions.
    string(APPEND DPLASMA_C_WFLAGS " -wd424,981,1419,1572,10237,11074,11076")
  ENDIF( CC_HAS_WD )
else(_match_icc)
  CHECK_C_COMPILER_FLAG( "-Wno-parentheses-equality" CC_HAS_PAR_EQUALITY )
  IF( CC_HAS_PAR_EQUALITY )
    string(APPEND DPLASMA_C_WFLAGS " -Wno-parentheses-equality")
  ENDIF( CC_HAS_PAR_EQUALITY )
endif(_match_icc)

# add gdb symbols in debug and relwithdebinfo
CHECK_C_COMPILER_FLAG( "-g3" CC_HAS_G3 )
IF( CC_HAS_G3 )
  string(APPEND CMAKE_C_FLAGS_DEBUG " -O0 -g3")
  string(APPEND CMAKE_C_FLAGS_RELWITHDEBINFO " -g3")
ELSE()
  string(APPEND CMAKE_C_FLAGS_DEBUG " -O0")
ENDIF( CC_HAS_G3)
# verbose compilation in debug
string(APPEND CMAKE_C_FLAGS_DEBUG " ${DPLASMA_C_WFLAGS}")
# remove asserts in release
string(APPEND CMAKE_C_FLAGS_RELEASE " -DNDEBUG")
string(APPEND CMAKE_C_FLAGS_RELWITHDEBINFO " ${DPLASMA_C_WFLAGS}")

#
# Fortran tricks: fortran main avoidance
#
IF (CMAKE_Fortran_COMPILER_WORKS)
  STRING(REGEX MATCH "Intel" _match_intel ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_intel)
    MESSAGE(STATUS "Add -nofor-main to the Fortran linker.")
    STRING(APPEND CMAKE_Fortran_LINK_EXECUTABLE " -nofor-main")
  ENDIF (_match_intel)

  STRING(REGEX MATCH "PGI$" _match_pgi ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_pgi)
    MESSAGE(STATUS "Add -Mnomain to the Fortran linker.")
    STRING(APPEND CMAKE_Fortran_LINK_EXECUTABLE " -Mnomain -Bstatic")
  ENDIF (_match_pgi)

  STRING(REGEX MATCH ".*xlf$" _match_xlf ${CMAKE_Fortran_COMPILER})
  IF (_match_xlf)
    MESSAGE(ERROR "Please use the thread-safe version of the xlf compiler (xlf_r)")
  ENDIF (_match_xlf)

#
# Even more Fortran tricks.
#
# Debug/Release FFLAGS depend on the compiler

  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
    MESSAGE(STATUS "Fortran adds libraries path ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
    FOREACH(ITEM IN ITEMS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
      list(APPEND EXTRA_LIBS "-L${ITEM}")
    ENDFOREACH(ITEM)
    MESSAGE(STATUS "Fortran adds libraries ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    list(APPEND EXTRA_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    # ifort
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
ENDIF (CMAKE_Fortran_COMPILER_WORKS)
