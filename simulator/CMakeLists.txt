FLEX_TARGET(dfllexer 
    dfllexer.l "${CMAKE_CURRENT_BINARY_DIR}/dfllexer.cpp" 
    COMPILE_FLAGS "-d"
)

BISON_TARGET(dflparser 
    dflparser.y "${CMAKE_CURRENT_BINARY_DIR}/dflparser.cpp" 
    DEFINES_FILE "${CMAKE_CURRENT_BINARY_DIR}/dflparser.h"
    COMPILE_FLAGS "-Wcounterexamples" 
)

ADD_FLEX_BISON_DEPENDENCY(dfllexer dflparser)

add_executable(sim 
    platform.h 
    platform.cpp 
    config.cpp
    location.h
    location.cpp
    main.cpp
    "${BISON_dflparser_OUTPUTS}"
    "${FLEX_dfllexer_OUTPUTS}"
    dflscanner.h
    parser.h
    parser.cpp
    parserextras.h
    progressbar.h
    parserextras.cpp 
    srccompiler.h 
    srccompiler.cpp 
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
    link_directories(${SuiteSparse_DIR}/lib)
    target_link_options(sim PUBLIC -static)
    target_link_libraries(
        sim 
        simlib
        Boost::filesystem
        ${SuiteSparse_DIR}/lib/libklu.a
        ${SuiteSparse_DIR}/lib/libamd.a
        ${SuiteSparse_DIR}/lib/libcamd.a
        ${SuiteSparse_DIR}/lib/libcolamd.a
        ${SuiteSparse_DIR}/lib/libccolamd.a
        ${SuiteSparse_DIR}/lib/libbtf.a
        ${SuiteSparse_DIR}/lib/libcholmod.a
        ${SuiteSparse_DIR}/lib/libklu_cholmod.a
        ${SuiteSparse_DIR}/lib/libsuitesparseconfig.a
        ${SuiteSparse_DIR}/lib/libopenblas.a
        m dl gomp wsock32
    )
else()
    # target_link_options(sim PUBLIC -static-libgcc -static-libstdc++)
    target_link_libraries(
        sim 
        simlib 
        dl 
        # -Wl,-Bstatic 
        m klu amd camd colamd ccolamd btf cholmod suitesparseconfig 
        boost_system boost_filesystem
        
    )
endif()

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
