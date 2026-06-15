cmake_minimum_required(VERSION 3.20)
project(hi_example LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)
find_package(Hyperiso REQUIRED)  # ~/.local/lib/cmake/Hyperiso

add_executable(core Core/base_core_example.cpp)
target_link_libraries(core PRIVATE Hyperiso::CoreLib)

add_executable(core_dep_pruner Core/dependency_pruner_example.cpp)
target_link_libraries(core_dep_pruner PRIVATE Hyperiso::CoreLib)

add_executable(wilson Wilson/base_wilson_example.cpp)
target_link_libraries(wilson PRIVATE Hyperiso::PhysicalModelLib)

add_executable(wilson_marty Wilson/wilson_with_marty_example.cpp)
target_link_libraries(wilson_marty PRIVATE Hyperiso::PhysicalModelLib)

add_executable(observable Observable/base_observable_example.cpp)
target_link_libraries(observable PRIVATE Hyperiso::BusinessLogicLib)

add_executable(chi2 Statistic/uncertainty_example.cpp)
target_link_libraries(chi2 PRIVATE Hyperiso::BusinessLogicLib)
