# Name of project
project(wignerSymbols)
set (wignerSymbols_VERSION_MAJOR 0)
set (wignerSymbols_VERSION_MINOR 1)
set (wignerSymbols_VERSION_RELEASE 1)

# CMake config
cmake_minimum_required(VERSION 2.8)
set (CMAKE_INSTALL_PREFIX /usr/)

# Compiler config
enable_language (Fortran)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -march=native")

# Included files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Source files
aux_source_directory(./src SRC_LIST)

# Build a shared library
add_library(${PROJECT_NAME} SHARED 
    ${SRC_LIST} 
    "./src/wignerSymbols-fortran.f"
    "./src/machine.for"
    "./src/fdump.f"
    "./src/j4save.f"
    "./src/xercnt.f"
    "./src/xerhlt.f"
    "./src/xermsg.f"
    "./src/xerprn.f"
    "./src/xersve.f"
    "./src/xgetua.f"
    "./src/wignerSymbols-fortran-c-binding.f90" )

SET_TARGET_PROPERTIES(${PROJECT_NAME}
    PROPERTIES
        VERSION ${wignerSymbols_VERSION_MAJOR}.${wignerSymbols_VERSION_MINOR}.${wignerSymbols_VERSION_RELEASE}
        SOVERSION ${wignerSymbols_VERSION_MAJOR}.${wignerSymbols_VERSION_MINOR}.${wignerSymbols_VERSION_RELEASE})

# Install directories
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
