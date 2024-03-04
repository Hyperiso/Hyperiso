cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


file(GLOB SOURCES "../Core/*.cpp" "../DataBase/*.cpp" "../logique_metier/*.cpp" "../Math/*.cpp" "../Physical_Model/*.cpp" "testParams.cpp")

add_executable(test ${SOURCES})

# Add include directories
target_include_directories(test 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
)

# Link any necessary libraries to the test executable
# For example:
# target_link_libraries(tests PRIVATE gtest)

# Optionally, add a custom target for running tests with CTest
# enable_testing()
# add_test(NAME AllTests COMMAND tests)

