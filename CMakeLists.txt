cmake_minimum_required(VERSION 3.18)

project(sim)

set(CMAKE_CXX_STANDARD 20)

enable_testing()

message(STATUS "CMake version ${CMAKE_VERSION}")
set(PROGRAM_NAME "vase")
set(PROGRAM_COPYRIGHT "2023-2024")
set(PROGRAM_VENDOR "EDA Lab FE Uni-Lj")

# Get version (versions in git are marked with tags of the form _<version number>)
# set(PROGRAM_VERSION "0.2")
execute_process(
    COMMAND git describe --tags --dirty --broken --match "_*"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
    RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (NOT GIT_DESCRIBE_ERROR_CODE)
    # Zero return code (OK)
    string(SUBSTRING "${GIT_DESCRIBE_VERSION}" 1 -1 PROGRAM_VERSION)
else()
    message(STATUS "Failed to retrieve version from git")
    set(PROGRAM_VERSION "unknown")
endif()
message(STATUS "Program version: ${PROGRAM_VERSION}")

# Platform
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
	# Windows
    set(SIM_PLATFORM "Windows")
    set(SIM_ARCH "x86_64")
	add_compile_definitions("SIMWINDOWS")
    set(PROGRAM_LIB_DIR "lib")
else()
	# Linux
    set(SIM_PLATFORM "Linux")
    set(SIM_ARCH "x86_64")
    set(PROGRAM_LIB_DIR "lib/${PROGRAM_NAME}")
endif()
message(STATUS "Platform: ${SIM_PLATFORM}")

# Compiler
if (DEFINED MSVC_VERSION)
	# Microsoft C, Windows
	message(STATUS "Using MSVC toolchain")
	add_compile_definitions("SIMMSC")
	# Floating point exceptions are disabled by default in MSVC
	# Need this for MSVC to correctly define __cplusplus, 
	# needed by Flex lexers and Bison parsers
	add_compile_options("/Zc:__cplusplus")
else()
	# GCC/MINGW
	message(STATUS "Using GCC/MINGW toolchain")
	add_compile_options("-fno-signaling-nans" "-fno-trapping-math")
endif()

set(STAGED_LIBRARY_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${PROGRAM_LIB_DIR}")
message(STATUS "Staged library directory: ${STAGED_LIBRARY_DIRECTORY}")

# Find Python (needed for tests)
find_package(Python 3 COMPONENTS Interpreter)

# Look for openvaf Verilog-A compiler
message(CHECK_START "Looking for openvaf Verilog-A compiler")
if (DEFINED OPENVAF_DIR)
  find_program(OPENVAF_COMPILER openvaf PATHS "${OPENVAF_DIR}" NO_DEFAULT_PATH)
else()
  find_program(OPENVAF_COMPILER openvaf)
endif()
if (OPENVAF_COMPILER STREQUAL "OPENVAF_COMPILER-NOTFOUND")
    message(FATAL_ERROR "  not found")
else()
    message(STATUS "  found ${OPENVAF_COMPILER}")    
endif()

#find_package(SWIG 4.0 COMPONENTS python)
#if(SWIG_FOUND)
#  include(UseSWIG)
#  message("SWIG found: ${SWIG_EXECUTABLE}")
#endif()

# find_package (Python3 COMPONENTS Interpreter Development)
# if(Python3_FOUND)
#   message("Python found: ${Python3_INCLUDE_DIRS}")
# endif()

find_package(Boost REQUIRED COMPONENTS filesystem)
include_directories("${Boost_INCLUDE_DIRS}")
link_directories("${Boost_LIBRARY_DIRS}")

#find_package(SuiteSparse CONFIG REQUIRED)
#include_directories("${SuiteSparse_INCLUDE_DIRS}")
include_directories("${SuiteSparse_DIR}/include")

find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)

# This is for circumventing a bug in mingw cmake 3.28 that causes it to incorrectly set two variables
if (${SIM_PLATFORM} STREQUAL "Windows")
	get_filename_component(tmp "${FLEX_EXECUTABLE}" DIRECTORY)
	set(FLEX_INCLUDE_DIRS "${tmp}/../include" CACHE STRING "Flex include dir" FORCE)
	set(FL_LIBRARY "${tmp}/../lib/libfl.a" CACHE STRING "Flex library" FORCE)
endif()

# Define SIMDEBUG preprocessor macro when build type is Debug
add_compile_definitions("$<$<CONFIG:Debug>:SIMDEBUG>")
# add_compile_definitions("$<$<CONFIG:Release>:SIMDEBUG>")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(lib)
add_subdirectory(simulator)
add_subdirectory(devices)
add_subdirectory(inc)
add_subdirectory(python)
add_subdirectory(test)

include("packaging.cmake")
