cmake_minimum_required(VERSION 3.20)
project(Hyperiso)
set(CMAKE_CXX_STANDARD 20)

set(SRC
    ./DataBase/lha_reader.cpp
    ./DataBase/lha_blocks.cpp
    ./DataBase/lha_elements.cpp
    ./Physical_Model/QCDParameters.cpp
    ./Core/Parameters.cpp
    ./Core/Logger.cpp 
    ./Core/MemoryManager.cpp
    ./Physical_Model/Wilson.cpp 
    ./Physical_Model/Wilson_parameters.cpp
    ./main_wilson.cpp
    ./Math/wilson_math.cpp
    ./Math/polylog.cpp 
    ./Math/others_functions.cpp 
    ./Math/exponential_integral.cpp 
    ./Math/bessel_function.cpp)

add_executable(hyperiso ${SRC})
