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


add_executable(demo2
    demo2.cpp
)

set_target_properties(demo2 PROPERTIES OUTPUT_NAME "demo2")
add_dependencies(demo2 simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(demo2 PUBLIC -static)
endif()

target_link_libraries(
    demo2
    simlib
    ${SIM_EXE_LIBRARIES}
)


add_executable(demo3
    demo3.cpp
)

set_target_properties(demo3 PROPERTIES OUTPUT_NAME "demo3")
add_dependencies(demo3 simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(demo3 PUBLIC -static)
endif()

target_link_libraries(
    demo3
    simlib
    ${SIM_EXE_LIBRARIES}
)


add_executable(demo4
    demo4.cpp
)

set_target_properties(demo4 PROPERTIES OUTPUT_NAME "demo4")
add_dependencies(demo4 simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(demo4 PUBLIC -static)
endif()

target_link_libraries(
    demo4
    simlib
    ${SIM_EXE_LIBRARIES}
)


add_executable(demo5
    demo5.cpp
)

set_target_properties(demo5 PROPERTIES OUTPUT_NAME "demo5")
add_dependencies(demo5 simlib)

if (${SIM_PLATFORM} STREQUAL "Windows")
    target_link_options(demo5 PUBLIC -static)
endif()

target_link_libraries(
    demo5
    simlib
    ${SIM_EXE_LIBRARIES}
)