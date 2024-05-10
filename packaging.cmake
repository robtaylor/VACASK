#
# Installation and packing
#

# This is the version of the product (<major>.<minor>)
set(CPACK_PACKAGE_NAME "VACASK")
set(CPACK_PACKAGE_VERSION "${PROGRAM_VERSION}")
set(CPACK_PACKAGE_VENDOR "${PROGRAM_VENDOR}")
set(CPACK_PACKAGE_CONTACT "arpad.buermen@fe.uni-lj.si")
    	
# Permissions
set(install_permissions_executable
	OWNER_READ OWNER_WRITE OWNER_EXECUTE
	GROUP_READ GROUP_EXECUTE
	WORLD_READ WORLD_EXECUTE
)

set(install_permisssions_file
	OWNER_READ OWNER_WRITE
	GROUP_READ
	WORLD_READ
)

set(install_permissions_directory ${install_permissions_executable})

# Set platform-dependent variables
if(${SIM_PLATFORM} STREQUAL "Windows")
    set(source_directory "src")
    set(docs_directory "doc")
else()
    set(source_directory "src/${PROGRAM_NAME}")
    set(docs_directory "share/doc/${PROGRAM_NAME}")
endif()

string(TOLOWER "${SIM_PLATFORM}" lowercase_platform)

# Executables
install(TARGETS sim
	DESTINATION "bin"
	PERMISSIONS ${install_permissions_executable}
)

install(FILES ${OPENVAF_COMPILER}
	DESTINATION "bin"
	PERMISSIONS ${install_permissions_executable}
)

# Device sources
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/devices" DESTINATION "${source_directory}"
	FILE_PERMISSIONS ${install_permissions_file}
	DIRECTORY_PERMISSIONS ${install_permissions_directory}
    PATTERN "CMakeLists.txt" EXCLUDE
)

# Library: compiled devices
foreach(item IN ITEMS ${OSDI_FILES})
    get_filename_component(tmpdirectory "${item}" DIRECTORY)
    install(FILES "${OSDI_DIR}/${item}" DESTINATION "${PROGRAM_LIB_DIR}/mod/${tmpdirectory}"
        PERMISSIONS ${install_permissions_file}
    )
endforeach()

# Library: simulator includes
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc" DESTINATION "${PROGRAM_LIB_DIR}"
	FILE_PERMISSIONS ${install_permissions_file}
	DIRECTORY_PERMISSIONS ${install_permissions_directory}
    FILES_MATCHING PATTERN "*.inc"
)

# Library: python files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python" DESTINATION "${PROGRAM_LIB_DIR}"
	FILE_PERMISSIONS ${install_permissions_file}
	DIRECTORY_PERMISSIONS ${install_permissions_directory}
    FILES_MATCHING PATTERN "*.py"
)

# License
install(FILES "LICENSE" DESTINATION "${docs_directory}"
    PERMISSIONS ${install_permissions_file}
)

# Demo files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/demo" DESTINATION "${docs_directory}"
	FILE_PERMISSIONS ${install_permissions_file}
	DIRECTORY_PERMISSIONS ${install_permissions_directory}
    PATTERN "CMakeLists.txt" EXCLUDE
)

# Test files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test" DESTINATION "${docs_directory}"
	FILE_PERMISSIONS ${install_permissions_file}
	DIRECTORY_PERMISSIONS ${install_permissions_directory}
    PATTERN "CMakeLists.txt" EXCLUDE
)


if(${SIM_PLATFORM} STREQUAL "Windows")
    set(CPACK_GENERATOR "ZIP")
else()
    set(CPACK_DEBIAN_PACKAGE_NAME "${PROGRAM_NAME}")
    set(CPACK_DEBIAN_DEBDROP_PACKAGE_NAME "ignore")

    set(CPACK_DEBIAN_PACKAGE_VERSION "${PROGRAM_VERSION}")

    if (SIM_ARCH STREQUAL "x86_64")
        set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    elseif (SPICEOPUS_PLATFORM STREQUAL "aarch64")
        set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "arm64")
    else()
        message(FATAL_ERROR "Unsupported architecture for building .deb package ${SIM_ARCH}")
    endif()
    message(STATUS "DEB package architecture: ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
    
    # Name of .deb file
    set(CPACK_DEBIAN_FILE_NAME "${CPACK_DEBIAN_PACKAGE_NAME}_${CPACK_DEBIAN_PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}.deb")
    message(STATUS "DEB package name: ${CPACK_DEBIAN_FILE_NAME}")

    # Generate shlib dependencies for deb package
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS "ON")

    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Arpad Buermen")
    set(CPACK_DEBIAN_PACKAGE_SECTION "electronics")
    set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Verilog-A Simulation Engine")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://fides.fe.uni-lj.si/vacask")

    # For Linux we build a tgz archive and a deb package
    set(CPACK_GENERATOR "TGZ" "DEB")
endif()

# Build package name
set(CPACK_PACKAGE_FILE_NAME "${PROGRAM_NAME}_${PROGRAM_VERSION}_${lowercase_platform}-${SIM_ARCH}")
message(STATUS "CPACK package name: ${CPACK_PACKAGE_FILE_NAME}")

include(CPack)
