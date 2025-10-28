add_executable(apidemo 
    apidemo.cpp
)

set_target_properties(apidemo PROPERTIES OUTPUT_NAME "apidemo")
target_include_directories(apidemo PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
target_link_libraries(apidemo ${CMAKE_CURRENT_BINARY_DIR}/../../lib/libsimlib.a)

add_dependencies(apidemo simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    link_directories(${SuiteSparse_DIR}/lib)
    target_link_options(apidemo PUBLIC -static)
    target_link_libraries(
        apidemo
        Boost::filesystem
        Boost::process
        Boost::system
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
        apidemo
        dl 
        # -Wl,-Bstatic 
        m klu amd camd colamd ccolamd btf cholmod suitesparseconfig 
        boost_system boost_filesystem boost_process
        
    )
endif()

