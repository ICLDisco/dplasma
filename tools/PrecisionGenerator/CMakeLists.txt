# Generate the Python3 precision generators if needed
# Note how this is done with an execute_process, during the configure
# rather than during the build: this is because GenDeps is invoked during
# configure, hence it needs to be created at this point
find_package(Python COMPONENTS Interpreter Development REQUIRED)
if(Python_VERSION_MAJOR GREATER 2)
  get_filename_component(PYTHON_EXE_DIR ${Python_EXECUTABLE} PATH)
  find_program(PYTHON_2TO3_EXECUTABLE
    NAMES 2to3
    HINTS ${PYTHON_EXE_DIR})
  if(NOT PYTHON_2TO3_EXECUTABLE)
    message(FATAL_ERROR "2to3 python utility not found. Use Python 2.7 or provide the 2to3 utility")
  endif()

  # If we are building in place and we need to convert the python tools, we
  # need to do the conversion in a temporary directory, to avoid overwriting the
  # original v2 tool
  if( DPLASMA_BUILD_INPLACE )
    set(PRECISIONGENERATOR_LOCATION ${CMAKE_BINARY_DIR}/tools/PrecisionGenerator/2to3)
    # we need to handle the case where the conversion was
    # already done once by removing the old files and start anew.
    if( EXISTS ${PRECISIONGENERATOR_LOCATION} )
      file(REMOVE_RECURSE ${PRECISIONGENERATOR_LOCATION})
    endif( EXISTS ${PRECISIONGENERATOR_LOCATION} )
  else( DPLASMA_BUILD_INPLACE )
    set(PRECISIONGENERATOR_LOCATION ${CMAKE_CURRENT_BINARY_DIR})
  endif( DPLASMA_BUILD_INPLACE )
  execute_process(OUTPUT_QUIET
    COMMAND ${PYTHON_2TO3_EXECUTABLE} -w -W -n -o ${PRECISIONGENERATOR_LOCATION} ${CMAKE_CURRENT_SOURCE_DIR})
else()
  set(PRECISIONGENERATOR_LOCATION ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set(PRECISIONGENERATOR_LOCATION ${PRECISIONGENERATOR_LOCATION} CACHE INTERNAL
  "Location of the PrecisionGenerator python tool to be created (should not be changed by users)")

install(FILES
        ${PRECISIONGENERATOR_LOCATION}/PrecisionDeps.py
        ${PRECISIONGENERATOR_LOCATION}/PrecisionGenerator.py
        ${PRECISIONGENERATOR_LOCATION}/subs.py
        ${PRECISIONGENERATOR_LOCATION}/Conversion.py
        DESTINATION ${CMAKE_INSTALL_LIBEXECDIR}/dplasma/PrecisionGenerator)
