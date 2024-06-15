cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(UserInterface)
add_executable(UserInterfaceLib main.cpp)

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target(UserInterfaceLib)
endif()
