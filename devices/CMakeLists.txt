set(device_sources
    "capacitor.va"
    "inductor.va"
    "resistor.va"
    "diode.va"
    "bsim3v3.va"
    "bsim4v8.va"
    "bsimbulk106.va"
)

# Prepare directory for staging the osdi library
set(STAGED_OSDI_DIRECTORY "${STAGED_LIBRARY_DIRECTORY}/mod")
file(MAKE_DIRECTORY "${STAGED_OSDI_DIRECTORY}")

# Set upt building and staging
set(osdi_built_files "")
set(osdi_staged_files "")
foreach(item IN ITEMS ${device_sources})
    get_filename_component(tmpnamebase "${item}" NAME_WLE)
    set(outputfile "${tmpnamebase}.osdi")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/build-${item}")
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${outputfile}"
        DEPENDS "${inputfile}"
        MAIN_DEPENDENCY "${inputfile}"
        COMMAND "${OPENVAF_COMPILER}" "-I${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${item}" -o "${outputfile}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${outputfile}" ".."
        COMMAND ${CMAKE_COMMAND} ARGS -E rm -f "${outputfile}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/build-${item}"
        COMMENT "Building osdi model from ${inputfile}"
    )
    list(APPEND osdi_built_files "${CMAKE_CURRENT_BINARY_DIR}/${outputfile}")
    add_custom_command(
        OUTPUT "${STAGED_OSDI_DIRECTORY}/${outputfile}"
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${outputfile}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_BINARY_DIR}/${outputfile}" "${STAGED_OSDI_DIRECTORY}/"
        COMMENT "Staging OSDI file ${outputfile}"
    )
    list(APPEND osdi_staged_files "${STAGED_OSDI_DIRECTORY}/${outputfile}")
endforeach()

add_custom_target(
    build_and_stage ALL 
    SOURCES "${device_sources}"
    DEPENDS "${osdi_staged_files}" "${osdi_built_files}"
)

SET(OSDI_FILES "${osdi_built_files}" CACHE STRING "List of OSDI files" FORCE)
