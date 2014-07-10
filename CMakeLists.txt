# Name of project
project(wignerSymbols)
set (wignerSymbols_VERSION_MAJOR 0)
set (wignerSymbols_VERSION_MINOR 0)
set (wignerSymbols_VERSION_RELEASE 1)

# CMake config
cmake_minimum_required(VERSION 2.8)
set (CMAKE_INSTALL_PREFIX /usr/)
enable_language (Fortran)

# Included files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Source files
aux_source_directory(./src SRC_LIST)

# Build a shared library
add_library(${PROJECT_NAME} SHARED ${SRC_LIST} "./src/wignerSymbols-fortran.f" "./src/wignerSymbols-fortran-c-binding.f90" )

SET_TARGET_PROPERTIES(${PROJECT_NAME}
    PROPERTIES
        VERSION ${wignerSymbols_VERSION_MAJOR}.${wignerSymbols_VERSION_MINOR}.${wignerSymbols_VERSION_RELEASE}
        SOVERSION ${wignerSymbols_VERSION_MAJOR}.${wignerSymbols_VERSION_MINOR}.${wignerSymbols_VERSION_RELEASE})

# Install directories
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
