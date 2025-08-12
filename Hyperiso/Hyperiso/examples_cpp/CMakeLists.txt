cmake_minimum_required(VERSION 3.20)
project(hi_example LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)
find_package(Hyperiso REQUIRED)  # utilise ~/.local/lib/cmake/Hyperiso

add_executable(core core_example.cpp)
target_link_libraries(core PRIVATE Hyperiso::CoreLib)

add_executable(wilson wilson_example.cpp)
target_link_libraries(wilson PRIVATE Hyperiso::PhysicalModelLib)

add_executable(observable observable_example.cpp)
target_link_libraries(observable PRIVATE Hyperiso::BusinessLogicLib)

add_executable(chi2 chi2_example.cpp)
target_link_libraries(chi2 PRIVATE Hyperiso::BusinessLogicLib)
