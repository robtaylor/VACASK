set(python_files
    rawfile.py
    runtest.py
    sg13g2tovc.py
    xschem2vc.py
    sg13g2xschem.tcl
    netlist_converter.py 
    netlist_converter/__init__.py 
    netlist_converter/converter.py 
    netlist_converter/dfl.py 
    netlist_converter/exc.py 
    netlist_converter/generators.py 
    netlist_converter/patterns.py
    netlist_converter/m_devices.py 
    netlist_converter/m_file.py 
    netlist_converter/m_inst_d.py 
    netlist_converter/m_inst_n.py 
    netlist_converter/m_inst_passive.py 
    netlist_converter/m_inst_q.py 
    netlist_converter/m_inst_x.py 
    netlist_converter/m_masters.py 
    netlist_converter/m_params.py
    netlist_converter/m_model.py 
    netlist_converter/m_output.py 
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
