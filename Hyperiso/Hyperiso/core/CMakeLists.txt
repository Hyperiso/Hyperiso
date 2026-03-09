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
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(FetchContent)

# define_property(GLOBAL PROPERTY SUPPORTED_OS BRIEF_DOC "Supported Operating Systems")
# set_property(GLOBAL PROPERTY SUPPORTED_OS "Linux")
# if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
#     message(FATAL_ERROR "This project is only supported on Linux.")
# endif()


if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MINUIT2)
    option(BUILD_WITH_MINUIT2 "Build project with Minuit2 support" OFF)
endif()

if (BUILD_WITH_2HDMC)
    add_compile_definitions(BUILD_WITH_2HDMC)
endif()

if (BUILD_WITH_SOFTSUSY)
    add_compile_definitions(BUILD_WITH_SOFTSUSY)
endif()

if (BUILD_WITH_SOFTSUSY)
    if(CMAKE_CXX_STANDARD LESS 14)
        message(FATAL_ERROR "SOFTSUSY requires at least C++14 support.")
    endif()

    find_program(MAKE_EXEC make PATHS /usr/local/bin /usr/bin /bin)
    if(NOT MAKE_EXEC)
        message(FATAL_ERROR "SOFTSUSY requires GNU Make. Install it before proceeding.")
    endif()

    find_program(FORTRAN_EXEC gfortran PATHS /usr/local/bin /usr/bin /bin)
    if(NOT FORTRAN_EXEC)
        message(FATAL_ERROR "SOFTSUSY requires a Fortran compiler (e.g., gfortran).")
    endif()

    find_program(LATEX_EXEC pdflatex PATHS /usr/local/bin /usr/bin /bin)
    if(NOT LATEX_EXEC)
        message(FATAL_ERROR "SOFTSUSY requires LaTeX. Install pdflatex and other LaTeX tools.")
    endif()
endif()

if (BUILD_WITH_2HDMC)
    message(STATUS "No specific dependency checks defined for 2HDMC in this script.")
endif()

if (BUILD_WITH_MARTY)
    find_program(CMAKE_EXEC cmake PATHS /usr/local/bin /usr/bin /bin)
    if(NOT CMAKE_EXEC)
        message(FATAL_ERROR "MARTY requires CMake. Install it before proceeding.")
    endif()

    find_program(LATEX_EXEC pdflatex PATHS /usr/local/bin /usr/bin /bin)
    if(NOT LATEX_EXEC)
        message(FATAL_ERROR "MARTY requires pdflatex for generating documentation.")
    endif()

    find_program(LUALATEX_EXEC lualatex PATHS /usr/local/bin /usr/bin /bin)
    if(NOT LUALATEX_EXEC)
        message(FATAL_ERROR "MARTY requires lualatex. Install it before proceeding.")
    endif()

    find_program(DVIPNG_EXEC dvipng PATHS /usr/local/bin /usr/bin /bin)
    if(NOT DVIPNG_EXEC)
        message(FATAL_ERROR "MARTY requires dvipng for documentation. Install dvipng.")
    endif()

    find_program(CONVERT_EXEC convert PATHS /usr/local/bin /usr/bin /bin)
    if(NOT CONVERT_EXEC)
        message(FATAL_ERROR "MARTY requires ImageMagick (convert tool). Install it before proceeding.")
    endif()
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

add_subdirectory(src/Common CommonLib)
add_subdirectory(src/Configured ConfiguredModule)
add_subdirectory(src/Math MathLib)
add_subdirectory(src/DataBase DataBaseLib)
add_subdirectory(src/ExternalIntegration ExternalIntegrationLib)
add_subdirectory(src/Core CoreLib)
add_subdirectory(src/ExternalIntegration/MartyInterface MartyInterfaceLib)
add_subdirectory(src/PhysicalModel PhysicalModelLib)
add_subdirectory(src/BusinessLogic BusinessLogicLib)
add_subdirectory(src/UserInterface UserInterfaceLib)

if (ENABLE_TESTS)
    # FetchContent_Declare(
    #     Catch2
    #     GIT_REPOSITORY https://github.com/catchorg/Catch2
    #     GIT_TAG v3.5.3
    #     GIT_SHALLOW TRUE)
    # FetchContent_MakeAvailable(Catch2)
    # list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    enable_testing()
    add_subdirectory(Test)
endif()

if (NOT DEFINED BUILD_WITH_PYTHON)
    option(BUILD_WITH_PYTHON "Build project with PYTHON support" OFF)
endif()

if (BUILD_WITH_PYTHON)
    add_subdirectory(src/Python pyhyperiso)
endif()


add_library(Hyperiso INTERFACE)

add_library(Hyperiso::Hyperiso ALIAS Hyperiso)

target_link_libraries(Hyperiso INTERFACE
    CommonLib
    ConfiguredModule
    MathLib
    DataBaseLib
    ExternalIntegrationLib
    CoreLib
    MartyInterfaceLib
    PhysicalModelLib
    BusinessLogicLib
    UserInterfaceLib
)

# Inclure src/ pendant le build (pour #include <hyperiso/...>)
target_include_directories(Hyperiso
    INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# --- Installation ---
include(GNUInstallDirs)

# Liste des modules à copier depuis src/ (temporaire tant que les publics sont dans src/)
set(HYPERISO_MODULES
    Common
    Configured
    Math
    DataBase
    ExternalIntegration
    Core
    ExternalIntegration/MartyInterface
    PhysicalModel
    BusinessLogic
    UserInterface
)

foreach(MOD IN LISTS HYPERISO_MODULES)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/${MOD}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/${MOD}
        FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "*.h"
            PATTERN "CMakeFiles" EXCLUDE
            PATTERN "*.cpp" EXCLUDE
    )
endforeach()

# Installe l’umbrella header
install(FILES ${PROJECT_SOURCE_DIR}/include/hyperiso.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso)

# # Installe la cible meta et les libs
# install(TARGETS
#     Hyperiso
#     CommonLib
#     ConfiguredModule
#     MathLib
#     DataBaseLib
#     ExternalIntegrationLib
#     CoreLib
#     MartyInterfaceLib
#     PhysicalModelLib
#     BusinessLogicLib
#     UserInterfaceLib
#     EXPORT HyperisoTargets
#     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )

install(TARGETS Hyperiso
    EXPORT HyperisoTargets
)

# Export pour find_package()
install(EXPORT HyperisoTargets
    NAMESPACE Hyperiso::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Hyperiso
)

include(CMakePackageConfigHelpers)

# Passe les options pour qu'elles soient substituées dans le .in
set(BUILD_WITH_SOFTSUSY ${BUILD_WITH_SOFTSUSY})
set(BUILD_WITH_2HDMC    ${BUILD_WITH_2HDMC})
set(BUILD_WITH_MARTY    ${BUILD_WITH_MARTY})

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/Cmake/HyperisoConfig.cmake.in   # <<< casse "Cmake"
    ${CMAKE_CURRENT_BINARY_DIR}/HyperisoConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Hyperiso
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/HyperisoConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/HyperisoConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/HyperisoConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Hyperiso
)