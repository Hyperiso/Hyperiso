cmake_minimum_required(VERSION 3.10)

project(TestRunner)


include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(test_runner test_runner.cpp)

add_test(NAME test_runner COMMAND test_runner)
