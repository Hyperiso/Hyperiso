cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(DataBase)

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()

add_library(DataBaseLib STATIC 
    lha_blocks.cpp
    lha_elements.cpp
    lha_reader.cpp)

set_target_properties(DataBaseLib PROPERTIES PUBLIC_HEADER "lha_elements.h;lha_blocks.h;lha_reader.h")

target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(DataBaseLib PUBLIC CommonLib)