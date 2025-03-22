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

if (BUILD_WITH_2HDMC)
    message(STATUS "wow")
    add_compile_definitions(BUILD_WITH_2HDMC)
endif()

if (BUILD_WITH_SOFTSUSY)
    add_definitions(-BUILD_WITH_SOFTSUSY)
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

if (NOT TARGET MartyInterfaceLib)
    add_subdirectory(../ExternalIntegration/MartyInterface MartyInterfaceLib)
endif()


add_library(PhysicalModelLib STATIC 
    domain/BWilsonRunningParameters.cpp
    domain/Wilson_parameters.cpp
    domain/Wilson.cpp
    domain/BWilson.cpp
    domain/WilsonGroup.cpp
    domain/WilsonManager.cpp
    adapters/ModelParamAdapter.cpp
    adapters/WilsonParamComposer.cpp
    BSM/epsilon_calculator.cpp
    BSM/susy_parameters.cpp
    BSM/thdm_parameters.cpp
    BSM/Wilson_SUSY.cpp
    BSM/Wilson_THDM.cpp)

set_target_properties(PhysicalModelLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(PhysicalModelLib PROPERTIES PUBLIC_HEADER "Interpolator.h;WilsonManager.h;Wilson.h;WilsonGroup.h;MartyWilson.h;Wilson_parameters.h;epsilon_calculator.h;susy_parameters.h;thdm_parameters;Wilson_susy.h;Wilson_THDM.h;WilsonManager.h")

target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/ports)
target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/domain)
target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/BSM)

target_link_libraries(PhysicalModelLib PUBLIC MartyInterfaceLib CoreLib CommonLib ConfiguredModule MathLib)

add_executable(test_wilson app/main_PM.cpp)
target_link_libraries(test_wilson PhysicalModelLib)

# add_executable(test_wilson main_phymod.cpp)
# target_link_libraries(test_wilson PhysicalModelLib)

# add_executable(test_martywilson main.cpp)
# target_link_libraries(test_martywilson PhysicalModelLib MartyInterfaceLib)
