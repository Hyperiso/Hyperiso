cmake_minimum_required(VERSION 3.20)
project(phyperiso)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

execute_process(
    COMMAND pybind11-config --cmakedir
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(PYBIND11_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
else()
    message(FATAL_ERROR "pybind11-config not found, please install pybind11.")
endif()


find_package(pybind11 CONFIG REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)


if (NOT TARGET CoreLib)
    add_subdirectory(../Core CoreLib)
endif()

if (NOT TARGET PhysicalModelLib)
    add_subdirectory(../PhysicalModel PhysicalModelLib)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/phyperiso)

add_subdirectory(src/core)
add_subdirectory(src/wilson)


# add_library(wilson_interface MODULE
#     src/wilson/binding_wilson.cpp
# )


