cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PhysicalModel)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()
if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured ConfiguredModule)
endif()
if (NOT TARGET MathLib)
    add_subdirectory(../Math MathLib)
endif()
if (NOT TARGET CoreLib)
    add_subdirectory(../Core CoreLib)
endif()

add_library(PhysicalModelLib STATIC 
    Wilson_parameters.cpp
    Wilson.cpp
    BSM/epsilon_calculator.cpp
    BSM/susy_parameters.cpp
    BSM/Wilson_susy.cpp
    BSM/Wilson_THDM.cpp)

set_target_properties(PhysicalModelLib PROPERTIES PUBLIC_HEADER "Wilson.h;Wilson_parameters.h;epsilon_calculator.h;susy_parameters.h;Wilson_susy.h;Wilson_THDM.h")

target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/BSM)

target_link_libraries(PhysicalModelLib PUBLIC CoreLib CommonLib ConfiguredModule MathLib)

add_executable(test_wilson Wilsonv2.cpp)
set_target_properties(test_wilson PROPERTIES PUBLIC_HEADER "Wilsonv2.h")
target_link_libraries(test_wilson PhysicalModelLib)
