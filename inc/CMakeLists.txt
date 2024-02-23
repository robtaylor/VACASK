set(inc_files
    builtins.inc
)

# Prepare directory for staging the osdi library
set(STAGED_INC_DIRECTORY "${STAGED_LIBRARY_DIRECTORY}/inc")
file(MAKE_DIRECTORY "${STAGED_INC_DIRECTORY}")

# Copy .py files to staging directory
set(inc_staged_files "")
foreach(item IN ITEMS ${inc_files})
    add_custom_command(
        OUTPUT "${STAGED_INC_DIRECTORY}/${item}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${item}"
        MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/${item}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${item}" "${STAGED_INC_DIRECTORY}"
        COMMENT "Staging simulator include file ${item}"
    )
    list(APPEND inc_staged_files "${STAGED_INC_DIRECTORY}/${item}")
endforeach()

add_custom_target(
    "stage_inc" ALL 
    SOURCES "${inc_files}"
    DEPENDS "${inc_staged_files}"
)
