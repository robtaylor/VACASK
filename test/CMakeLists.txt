# Testing - requires a working installation of python3 and openvaf

# Directory where the python script that runs a test can be found
# Use staging directory to keep source clean
set(PYTHON_LIB "${STAGED_LIBRARY_DIRECTORY}/python")
message("Python library path for testing: ${PYTHON_LIB}")

# macro that sets up a test
macro(prepare_netlist_test simfile)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${simfile}")
    add_test(NAME "${simfile}"
        COMMAND $<TARGET_FILE:sim> -dp "${CMAKE_CURRENT_SOURCE_DIR}/${simfile}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${simfile}"
    )
    
    # message("Adding netlist test ${simfile}")
    
    # PYTHONPATH points to the binary test directory
    # Each test has a subdirectory there
    # Dependencies (.py files like diocalc.py) will be copied there
    set_tests_properties(
        "${simfile}" PROPERTIES 
        ENVIRONMENT 
        "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR};SIM_TEST=yes"
        TIMEOUT 5
    )
endmacro()

# Create directory for manually running tests from vscode
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/manual_run")

set(TEST_DEPENDENCIES
    diocalc.py
)

foreach(item IN ITEMS ${TEST_DEPENDENCIES})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${item} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${item} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/manual_run")
endforeach()

prepare_netlist_test(test_op.sim)
prepare_netlist_test(test_opnl.sim)
prepare_netlist_test(test_opnodeset.sim)
prepare_netlist_test(test_opstore.sim)
prepare_netlist_test(test_opsweep.sim)
prepare_netlist_test(test_sweep.sim)
prepare_netlist_test(test_sweeppar.sim)
prepare_netlist_test(test_sweepvar.sim)
prepare_netlist_test(test_sweepvec.sim)
prepare_netlist_test(test_dcinc.sim)
prepare_netlist_test(test_dctf.sim)
prepare_netlist_test(test_ac.sim)
prepare_netlist_test(test_actf.sim)
prepare_netlist_test(test_tran.sim)
prepare_netlist_test(test_tranic.sim)
prepare_netlist_test(test_noise.sim)
prepare_netlist_test(test_build.sim)
prepare_netlist_test(test_visrc.sim)
prepare_netlist_test(test_viwfm.sim)
prepare_netlist_test(test_ctlsrc.sim)
prepare_netlist_test(test_resistor.sim)
prepare_netlist_test(test_capacitor.sim)
prepare_netlist_test(test_inductor.sim)
prepare_netlist_test(test_mutual.sim)
prepare_netlist_test(test_diode.sim)
prepare_netlist_test(test_inverter.sim)
