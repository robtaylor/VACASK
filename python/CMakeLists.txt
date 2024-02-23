set(python_files
    rawfile.py
    runtest.py
)

# Prepare directory for staging the osdi library
set(STAGED_PYTHON_DIRECTORY "${STAGED_LIBRARY_DIRECTORY}/python")
file(MAKE_DIRECTORY "${STAGED_PYTHON_DIRECTORY}")

# Copy .py files to staging directory
set(python_staged_files "")
foreach(item IN ITEMS ${python_files})
    add_custom_command(
        OUTPUT "${STAGED_PYTHON_DIRECTORY}/${item}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${item}"
        MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/${item}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${item}" "${STAGED_PYTHON_DIRECTORY}"
        COMMENT "Staging Python file ${item}"
    )
    list(APPEND python_staged_files "${STAGED_PYTHON_DIRECTORY}/${item}")
endforeach()

add_custom_target(
    "stage_scripts" ALL 
    SOURCES "${python_files}"
    DEPENDS "${python_staged_files}"
)
