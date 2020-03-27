# Generate the Python3 precision generators if needed
# Note how this is done with an execute_process, during the configure
# rather than during the build: this is because GenDeps is invoked during
# configure, hence it needs to be created at this point
find_package(Python REQUIRED COMPONENTS Interpreter Development)
if(Python_VERSION_MAJOR GREATER 2)
  get_filename_component(PYTHON_EXE_DIR ${Python_EXECUTABLE} PATH)
  find_program(PYTHON_2TO3_EXECUTABLE
    NAMES 2to3
    HINTS ${PYTHON_EXE_DIR})
  if(NOT PYTHON_2TO3_EXECUTABLE)
    message(FATAL_ERROR "2to3 python utility not found. Use Python 2.7 or provide the 2to3 utility")
  endif()

  execute_process(OUTPUT_QUIET
      COMMAND ${PYTHON_2TO3_EXECUTABLE} -w -W -n -o ${CMAKE_BINARY_DIR}/tools/PrecisionGenerator ${CMAKE_SOURCE_DIR}/tools/PrecisionGenerator)

    install(FILES
            ${CMAKE_BINARY_DIR}/tools/PrecisionGenerator/PrecisionDeps.py
            ${CMAKE_BINARY_DIR}/tools/PrecisionGenerator/PrecisionGenerator.py
            ${CMAKE_BINARY_DIR}/tools/PrecisionGenerator/subs.py
            ${CMAKE_BINARY_DIR}/tools/PrecisionGenerator/Conversion.py
            DESTINATION ${CMAKE_INSTALL_LIBEXECDIR}/dplasma/PrecisionGenerator)
else()
    install(FILES
            ${CMAKE_SOURCE_DIR}/tools/PrecisionGenerator/PrecisionDeps.py
            ${CMAKE_SOURCE_DIR}/tools/PrecisionGenerator/PrecisionGenerator.py
            ${CMAKE_SOURCE_DIR}/tools/PrecisionGenerator/subs.py
            ${CMAKE_SOURCE_DIR}/tools/PrecisionGenerator/Conversion.py
            DESTINATION ${CMAKE_INSTALL_LIBEXECDIR}/dplasma/PrecisionGenerator)
endif()

