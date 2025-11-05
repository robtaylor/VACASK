add_executable(apidemo 
    apidemo.cpp
)

set_target_properties(apidemo PROPERTIES OUTPUT_NAME "apidemo")
add_dependencies(apidemo simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(apidemo PUBLIC -static)
endif()

target_link_libraries(
    apidemo
    simlib
    ${SIM_EXE_LIBRARIES}
)
