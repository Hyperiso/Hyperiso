cmake_minimum_required(VERSION 3.20)
project(Hyperiso VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Cmake/")
include(Docs)
include(Warnings)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

option(ENABLE_TESTS "Enable building of tests" OFF)

option(ENABLE_WARNINGS "Enable warnings" ON)
option(ENABLE_WARNINGS_AS_ERRORS "Enable warnings as errors" OFF)

option(ENABLE_ADRESS_SANITIZER "Enable adress sanitizer" ON)
option(ENABLE_UNDEFINED_SANITIZER "Enable undefined sanitizer" ON)

if (ENABLE_ADRESS_SANITIZER OR ENABLE_UNDEFINED_SANITIZER)
    include(Sanitizers)
    add_sanitizer_flags()
endif()

add_subdirectory(Common CommonLib)
add_subdirectory(Configured ConfiguredModule)
add_subdirectory(Math MathLib)
add_subdirectory(DataBase DataBaseLib)
add_subdirectory(Core CoreLib)
add_subdirectory(PhysicalModel PhysicalModelLib)
add_subdirectory(BusinessLogic BusinessLogicLib)
add_subdirectory(UserInterface UserInterfaceLib)

if (ENABLE_TESTS)
    add_subdirectory(Test)
endif()