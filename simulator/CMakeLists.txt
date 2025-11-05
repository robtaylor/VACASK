add_executable(sim 
    platform.h 
    platform.cpp 
    config.cpp
    main.cpp
    progressbar.h
    progressbar.cpp 
    cmd.h 
    cmd.cpp 
)
set_target_properties(sim PROPERTIES OUTPUT_NAME "${PROGRAM_NAME}")
target_include_directories(sim PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "${FLEX_INCLUDE_DIRS}")

# Program information
target_compile_definitions(sim PUBLIC PROGRAM_NAME=${PROGRAM_NAME})
target_compile_definitions(sim PUBLIC PROGRAM_VERSION=${PROGRAM_VERSION})
target_compile_definitions(sim PUBLIC PROGRAM_COPYRIGHT=${PROGRAM_COPYRIGHT})
target_compile_definitions(sim PUBLIC PROGRAM_HOMEPAGE="${PROGRAM_HOMEPAGE}")

add_dependencies(sim simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(sim PUBLIC -static)
endif()

target_link_libraries(
    sim 
    simlib 
    ${SIM_EXE_LIBRARIES}
)

# Stage openvaf compiler
get_filename_component(openvaf_executable_name "${OPENVAF_COMPILER}" NAME)
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${openvaf_executable_name}"
    DEPENDS "${OPENVAF_COMPILER}"
    COMMAND ${CMAKE_COMMAND} ARGS -E copy "${OPENVAF_COMPILER}" "${CMAKE_CURRENT_BINARY_DIR}/${openvaf_executable_name}"
    COMMENT "Staging openvaf compiler as: ${openvaf_executable_name}"
)
add_custom_target(
    "stage_openvaf" ALL 
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${openvaf_executable_name}"
)
