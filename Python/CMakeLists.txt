cmake_minimum_required(VERSION 3.20)
project(WilsonInterfaceProject)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(pybind11 REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)

if (NOT TARGET PhysicalModelLib)
    add_subdirectory(../PhysicalModel PhysicalModelLib)
endif()

add_library(wilson_interface MODULE
    src/wilson/binding_wilson.cpp
    src/wilson/binding_manager.cpp
    src/wilson/binding_utilities.cpp
)

target_link_libraries(wilson_interface PRIVATE pybind11::module PhysicalModelLib)

set_target_properties(wilson_interface PROPERTIES
    PREFIX ""
    SUFFIX ".so"
)
