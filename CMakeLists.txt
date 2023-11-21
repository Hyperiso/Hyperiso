cmake_minimum_required(VERSION 3.20)
project(Hyperiso)
set(CMAKE_CXX_STANDARD 17)

set(SRC
    ./DataBase/lha_elements.cpp
    ./DataBase/lha_reader.cpp)

add_executable(hyperiso ${SRC})