cmake_minimum_required(VERSION 3.20)
project(phyperiso LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# pybind11 via pybind11-config (OK tant que le package pybind11 est installé via pip)
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

if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured ConfiguredModule)
endif()

if (NOT TARGET CoreLib)
    add_subdirectory(../Core CoreLib)
endif()

if (NOT TARGET PhysicalModelLib)
    add_subdirectory(../PhysicalModel PhysicalModelLib)
endif()

if (NOT TARGET BusinessLogicLib)
    add_subdirectory(../BusinessLogic BusinessLogicLib)
endif()

if (NOT TARGET StatisticLib)
    add_subdirectory(../Statistic StatisticLib)
endif()

add_library(pyhyperiso MODULE
    src/pyhyperiso.cpp
    src/common/binding_common.cpp
    src/math/binding_math.cpp
    src/core/binding_core.cpp
    src/wilson/binding_wilson.cpp
    src/observable/binding_observable.cpp
    src/statistic/binding_statistic.cpp
)

target_link_libraries(pyhyperiso PRIVATE
    pybind11::module
    CoreLib
    PhysicalModelLib
    BusinessLogicLib
    StatisticLib
)

set_target_properties(pyhyperiso PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)

# ⚠️ IMPORTANT : Installation dans le wheel
# Destination RELATIVE à la racine de site-packages
# => le .so sera installé dans pyhyperiso/phyperiso/ dans le wheel
install(TARGETS pyhyperiso
        LIBRARY DESTINATION "pyhyperiso/phyperiso")
