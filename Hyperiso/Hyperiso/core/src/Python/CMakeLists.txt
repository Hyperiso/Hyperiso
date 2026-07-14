cmake_minimum_required(VERSION 3.20)
if(DEFINED SKBUILD_PROJECT_VERSION)
    project(phyperiso VERSION "${SKBUILD_PROJECT_VERSION}" LANGUAGES CXX)
else()
    project(phyperiso VERSION 1.0.0 LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
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

pybind11_add_module(pyhyperiso MODULE
    src/pyhyperiso.cpp
    src/common/binding_common.cpp
    src/math/binding_math.cpp
    src/core/binding_core.cpp
    src/wilson/binding_wilson.cpp
    src/observable/binding_observable.cpp
    src/statistic/binding_statistic.cpp
)

target_link_libraries(pyhyperiso PRIVATE
    CoreLib
    PhysicalModelLib
    BusinessLogicLib
    StatisticLib
)


install(TARGETS pyhyperiso
        LIBRARY DESTINATION "pyhyperiso/phyperiso")
