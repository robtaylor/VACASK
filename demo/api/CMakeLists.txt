add_executable(demo1
    demo1.cpp
)

set_target_properties(demo1 PROPERTIES OUTPUT_NAME "demo1")
add_dependencies(demo1 simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(demo1 PUBLIC -static)
endif()

target_link_libraries(
    demo1
    simlib
    ${SIM_EXE_LIBRARIES}
)
