# Name of project
project(wignerSymbols)
enable_testing()

set (wignerSymbols_VERSION_MAJOR 0)
set (wignerSymbols_VERSION_MINOR 2)
set (wignerSymbols_VERSION_RELEASE 0)

# CMake config
cmake_minimum_required(VERSION 2.8)
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX /usr)
endif()

# Compiler config
enable_language (Fortran)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wpedantic -Werror -march=native")

# Coverage report
option(WITH_COVERAGE "Generate code coverage report" OFF)
if(WITH_COVERAGE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  link_libraries(gcov)
endif()

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

add_subdirectory(tests)

# Install directories
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
