set(python_files
    rawfile.py
    runtest.py
    sg13g2tovc.py
    ng2vc.py 
    ng2vclib/__init__.py 
    ng2vclib/converter.py 
    ng2vclib/dfl.py 
    ng2vclib/exc.py 
    ng2vclib/generators.py 
    ng2vclib/m_devices.py 
    ng2vclib/m_file.py 
    ng2vclib/m_inst_d.py 
    ng2vclib/m_inst_n.py 
    ng2vclib/m_inst_passive.py 
    ng2vclib/m_inst_q.py 
    ng2vclib/m_masters.py 
    ng2vclib/m_model.py 
    ng2vclib/m_output.py 
)

# Prepare directory for staging the osdi library
set(STAGED_PYTHON_DIRECTORY "${STAGED_LIBRARY_DIRECTORY}/python")
file(MAKE_DIRECTORY "${STAGED_PYTHON_DIRECTORY}")

# Copy .py files to staging directory
set(python_staged_files "")
foreach(item IN ITEMS ${python_files})
    # Get directory name
    get_filename_component(tmpdirectory "${item}" DIRECTORY)
    if (NOT "${tmpdirectory}" STREQUAL "")
        file(MAKE_DIRECTORY "${STAGED_PYTHON_DIRECTORY}/${tmpdirectory}")
    endif()
    add_custom_command(
        OUTPUT "${STAGED_PYTHON_DIRECTORY}/${item}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${item}"
        MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/${item}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${item}" "${STAGED_PYTHON_DIRECTORY}/${tmpdirectory}"
        COMMENT "Staging Python file ${item}"
    )
    list(APPEND python_staged_files "${STAGED_PYTHON_DIRECTORY}/${item}")
endforeach()

add_custom_target(
    "stage_scripts" ALL 
    SOURCES "${python_files}"
    DEPENDS "${python_staged_files}"
)
