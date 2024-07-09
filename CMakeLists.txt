cmake_minimum_required(VERSION 3.20)
project(Hyperiso VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Cmake/")
include(Docs)
include(Warnings)
include(Tools)

include(FetchContent)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()

option(ENABLE_TESTS "Enable building of tests" OFF)

option(ENABLE_WARNINGS "Enable warnings" OFF)
option(ENABLE_WARNINGS_AS_ERRORS "Enable warnings as errors" OFF)

option(ENABLE_ADRESS_SANITIZER "Enable adress sanitizer" OFF)
option(ENABLE_UNDEFINED_SANITIZER "Enable undefined sanitizer" OFF)

option(ENABLE_CLANG_TIDY "Enable to add clang tidy." OFF)

option(ENABLE_CLANG_FORMAT "Enable to add clang-format." OFF)
option(ENABLE_CMAKE_FORMAT "Enable to add cmake-format." OFF)


if (ENABLE_ADRESS_SANITIZER OR ENABLE_UNDEFINED_SANITIZER)
    include(Sanitizers)
    add_sanitizer_flags()
endif()

add_cmake_format_target()
add_clang_format_target()

add_subdirectory(Common CommonLib)
add_subdirectory(Configured ConfiguredModule)
add_subdirectory(Math MathLib)
add_subdirectory(DataBase DataBaseLib)
add_subdirectory(ExternalIntegration ExternalIntegrationLib)
add_subdirectory(Core CoreLib)
add_subdirectory(PhysicalModel PhysicalModelLib)
add_subdirectory(BusinessLogic BusinessLogicLib)
add_subdirectory(UserInterface UserInterfaceLib)

if (ENABLE_TESTS)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2
        GIT_TAG v3.5.3
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    enable_testing()
    add_subdirectory(Test)
endif()