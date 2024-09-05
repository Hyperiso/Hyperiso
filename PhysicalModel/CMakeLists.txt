cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 23)
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
    Wilsonv2.cpp
    WilsonManager.cpp
    BSM/epsilon_calculator.cpp
    BSM/susy_parameters.cpp
    BSM/thdm_parameters.cpp
    BSM/Wilson_susyv2.cpp
    BSM/Wilson_THDMv2.cpp)

set_target_properties(PhysicalModelLib PROPERTIES PUBLIC_HEADER "Wilsonv2.h;Wilson_parameters.h;epsilon_calculator.h;susy_parameters.h;thdm_parameters;Wilson_susyv2.h;Wilson_THDMv2.h;WilsonManager.h")

target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/BSM)

target_link_libraries(PhysicalModelLib PUBLIC CoreLib CommonLib ConfiguredModule MathLib)

add_executable(test_wilson main_phymod.cpp)
# set_target_properties(test_wilson PROPERTIES PUBLIC_HEADER "Wilsonv2.h")
target_link_libraries(test_wilson PhysicalModelLib)
