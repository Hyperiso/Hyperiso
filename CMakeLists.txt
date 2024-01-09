cmake_minimum_required(VERSION 3.20)
project(Hyperiso)
set(CMAKE_CXX_STANDARD 20)

set(SRC
    ./Physical_Model/QCDParameters.cpp
    ./main.cpp)

add_executable(hyperiso ${SRC})
