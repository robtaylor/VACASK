cmake_minimum_required(VERSION 3.18)

project(sim)

set(CMAKE_CXX_STANDARD 20)
if(POLICY CMP0177)
    cmake_policy(SET CMP0177 NEW)
endif()

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

enable_testing()

message(STATUS "CMake version ${CMAKE_VERSION}")
set(PROGRAM_NAME "vacask")
set(PROGRAM_COPYRIGHT "2023-2025")
set(PROGRAM_VENDOR "EDA Lab FE Uni-Lj")
set(PROGRAM_HOMEPAGE "https://codeberg.org/arpadbuermen/VACASK")

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
	add_compile_options("-fno-signaling-nans" "-fno-trapping-math" "-fcoroutines")
endif()

set(STAGED_LIBRARY_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${PROGRAM_LIB_DIR}")
message(STATUS "Staged library directory: ${STAGED_LIBRARY_DIRECTORY}")

# Find Python (needed for tests)
find_package(Python 3 COMPONENTS Interpreter)

# Look for openvaf Verilog-A compiler
message(CHECK_START "Looking for OpenVAF-reloaded Verilog-A compiler")
if (DEFINED OPENVAF_DIR)
  find_program(OPENVAF_COMPILER openvaf-r PATHS "${OPENVAF_DIR}" NO_DEFAULT_PATH)
else()
  find_program(OPENVAF_COMPILER openvaf-r)
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

if (${SIM_PLATFORM} STREQUAL "Windows")
    find_package(Boost REQUIRED COMPONENTS filesystem process)
    include_directories("${Boost_INCLUDE_DIRS}")
    link_directories("${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost include: ${Boost_INCLUDE_DIRS}")
    include_directories("${TOMLPP_DIR}/include")
    message(STATUS "TOMLPP include: ${TOMLPP_DIR}/include")
else()
    # For 1.83
    # set(CMAKE_PREFIX_PATH "/usr/lib/x86_64-linux-gnu/cmake/Boost-1.88.0")
    # find_package(Boost CONFIG REQUIRED filesystem)

    # For 1.88 installed manually by user
    set(Boost_NO_SYSTEM_PATHS TRUE)
    find_package(Boost 1.88 REQUIRED COMPONENTS filesystem process system)
    include_directories("${Boost_INCLUDE_DIRS}")
    set(Boost_EXTRA_LINK_DIR "${Boost_INCLUDE_DIRS}/stage/lib")
    link_directories(${Boost_EXTRA_LINK_DIR})
    message(STATUS "Boost include: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost extra link directory: ${Boost_EXTRA_LINK_DIR}")

    # For 1.88 system-installed boost
    # find_package(Boost 1.88 REQUIRED COMPONENTS filesystem process)
    # include_directories("${Boost_INCLUDE_DIRS}")
    # link_directories("${Boost_LIBRARY_DIRS}")
    # message(STATUS "Boost include: ${Boost_INCLUDE_DIRS} lib ${Boost_EXTRA_LINK_DIR}")
endif()

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

#
# Custom build types
# 

# Profiling debug version -pg -fprofile-arcs -ftest-coverage --coverage
# Also does line coverage
# Function profile
#   gprof <binary> gmon.out
# Line-by-line profile
#   gprof -l <binary> gmon.out
set(CMAKE_CXX_FLAGS_DEBUGPROFILE "${CMAKE_CXX_FLAGS_DEBUG} -pg -fprofile-arcs -ftest-coverage --coverage" CACHE STRING
    "Flags used by the C++ compiler during debug profile builds." FORCE)
set(CMAKE_C_FLAGS_DEBUGPROFILE "${CMAKE_C_FLAGS_DEBUG} -pg -fprofile-arcs -ftest-coverage --coverage" CACHE STRING
    "Flags used by the C compiler during debug profile builds." FORCE)
set(CMAKE_EXE_LINKER_FLAGS_DEBUGPROFILE "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg" CACHE STRING
    "Flags used for linking binaries during debug profile builds." FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_DEBUGPROFILE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -pg" CACHE STRING
    "Flags used by the shared libraries linker during debug profile builds." FORCE)

# Profiling release version
set(CMAKE_CXX_FLAGS_RELEASEPROFILE "${CMAKE_CXX_FLAGS_RELEASE} -pg -fprofile-arcs -ftest-coverage --coverage" CACHE STRING
"Flags used by the C++ compiler during release profile builds." FORCE)
set(CMAKE_C_FLAGS_RELEASEPROFILE "${CMAKE_C_FLAGS_RELEASE} -pg -fprofile-arcs -ftest-coverage --coverage" CACHE STRING
"Flags used by the C compiler during release profile builds." FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELEASEPROFILE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -pg" CACHE STRING
"Flags used for linking binaries during release profile builds." FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELEASEPROFILE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -pg" CACHE STRING
"Flags used by the shared libraries linker during release profile builds." FORCE)

# Address sanitizer debug version
set(CMAKE_CXX_FLAGS_ASANDEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address" CACHE STRING
"Flags used by the C++ compiler during address sanitizer debug builds." FORCE)
set(CMAKE_C_FLAGS_ASANDEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address" CACHE STRING
"Flags used by the C compiler during address sanitizer debug builds." FORCE)
set(CMAKE_EXE_LINKER_FLAGS_ASANDEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -static-libasan" CACHE STRING
"Flags used for linking binaries during address sanitizer debug builds." FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_ASANDEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -static-libasan" CACHE STRING
"Flags used by the shared libraries linker during address sanitizer debug builds." FORCE)

MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_DEBUGPROFILE
    CMAKE_C_FLAGS_DEBUGPROFILE
    CMAKE_EXE_LINKER_FLAGS_DEBUGPROFILE
    CMAKE_SHARED_LINKER_FLAGS_DEBUGPROFILE

    CMAKE_CXX_FLAGS_RELEASEPROFILE
    CMAKE_C_FLAGS_RELEASEPROFILE
    CMAKE_EXE_LINKER_FLAGS_RELEASEPROFILE
    CMAKE_SHARED_LINKER_FLAGS_RELEASEPROFILE

    CMAKE_CXX_FLAGS_ASANDEBUG
    CMAKE_C_FLAGS_ASANDEBUG
    CMAKE_EXE_LINKER_FLAGS_ASANDEBUG
    CMAKE_SHARED_LINKER_FLAGS_ASANDEBUG
)

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Debug
       CACHE STRING "Choose the type of build : None Debug Release DebugProfile ReleaseProfile AsanDebug."
       FORCE
    )
ENDIF(NOT CMAKE_BUILD_TYPE)
message("* Current build type is : ${CMAKE_BUILD_TYPE}")

#
# End of custom build types
# 


# Define SIMDEBUG preprocessor macro when build type is Debug
add_compile_definitions("$<$<CONFIG:Debug>:SIMDEBUG>")
add_compile_definitions("$<$<CONFIG:DebugProfile>:SIMDEBUG>")
add_compile_definitions("$<$<CONFIG:ReleaseProfile>:SIMPROFILE>")
add_compile_definitions("$<$<CONFIG:DebugProfile>:SIMPROFILE>")
add_compile_definitions("$<$<CONFIG:AsanDebug>:SIMDEBUG>")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(lib)
add_subdirectory(simulator)
add_subdirectory(devices)
add_subdirectory(inc)
add_subdirectory(python)
add_subdirectory(test)

add_subdirectory(demo/api)

include("packaging.cmake")
