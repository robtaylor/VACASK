set(device_sources
    "capacitor.va"
    "inductor.va"
    "resistor.va"
    "diode.va"
    "bsim3v3.va"
    "bsim4v8.va"
    "bsimbulk106.va"
    "spice/resistor.va"
    "spice/capacitor.va"
    "spice/inductor.va"
    "spice/diode.va"
    "spice/jfet.va"
    "spice/bjt.va"
)

# Prepare directory for staging the osdi library
set(STAGED_OSDI_DIRECTORY "${STAGED_LIBRARY_DIRECTORY}/mod")
file(MAKE_DIRECTORY "${STAGED_OSDI_DIRECTORY}")

# Set upt building and staging
set(osdi_built_files "")
set(osdi_staged_files "")
set(osdi_subdirs "")
set(osdi_built_files_rel)
foreach(item IN ITEMS ${device_sources})
    set(inputfile "${CMAKE_CURRENT_SOURCE_DIR}/${item}")
    get_filename_component(tmpnamebase "${item}" NAME_WLE)
    get_filename_component(tmpdirectory "${item}" DIRECTORY)
    if ("${tmpdirectory}" EQUAL "")
        set(tmpdirectory ".")
        set(builddir "${CMAKE_CURRENT_BINARY_DIR}/build-${tmpnamebase}")
    else()
        set(builddir "${CMAKE_CURRENT_BINARY_DIR}/${tmpdirectory}/build-${tmpnamebase}")
    endif()
    if (NOT "${tmpdirectory}" EQUAL "." AND NOT "${tmpdirectory}" IN_LIST osdi_subdirs)
        list(APPEND osdi_subdirs "${tmpdirectory}")
        file(MAKE_DIRECTORY "${STAGED_OSDI_DIRECTORY}/${tmpdirectory}")
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${tmpdirectory}")
    endif()
    message(STATUS ${builddir})
    set(outputfile "${tmpnamebase}.osdi")
    file(MAKE_DIRECTORY "${builddir}")
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${tmpdirectory}/${outputfile}"
        DEPENDS "${inputfile}"
        MAIN_DEPENDENCY "${inputfile}"
        COMMAND "${OPENVAF_COMPILER}" "--allow" "variant_const_simparam" "-I${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${item}" -o "${outputfile}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${outputfile}" ".."
        COMMAND ${CMAKE_COMMAND} ARGS -E rm -f "${outputfile}"
        WORKING_DIRECTORY "${builddir}"
        COMMENT "Building osdi model from ${inputfile}"
    )
    list(APPEND osdi_built_files "${CMAKE_CURRENT_BINARY_DIR}/${tmpdirectory}/${outputfile}")
    list(APPEND osdi_built_files_rel "${tmpdirectory}/${outputfile}")
    add_custom_command(
        OUTPUT "${STAGED_OSDI_DIRECTORY}/${tmpdirectory}/${outputfile}"
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${tmpdirectory}/${outputfile}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_BINARY_DIR}/${tmpdirectory}/${outputfile}" "${STAGED_OSDI_DIRECTORY}/${tmpdirectory}/"
        COMMENT "Staging OSDI file ${outputfile}"
    )
    list(APPEND osdi_staged_files "${STAGED_OSDI_DIRECTORY}/${tmpdirectory}/${outputfile}")
endforeach()

add_custom_target(
    build_and_stage_devices ALL 
    SOURCES "${device_sources}"
    DEPENDS "${osdi_staged_files}" "${osdi_built_files}"
)

SET(OSDI_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE STRING "Directory where OSDI files are atored" FORCE)
SET(OSDI_FILES "${osdi_built_files_rel}" CACHE STRING "List of OSDI files relative to OSDI_DIR" FORCE)
