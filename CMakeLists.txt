cmake_minimum_required(VERSION 3.20)
project(Hyperiso)
set(CMAKE_CXX_STANDARD 20)

set(SRC
    ./DataBase/lha_reader.cpp
    ./DataBase/lha_blocks.cpp
    ./DataBase/lha_elements.cpp)
    ./Physical_Model/QCDParameters.cpp
    ./main.cpp)

add_executable(hyperiso ${SRC})
