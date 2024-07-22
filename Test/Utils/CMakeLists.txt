cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Utils)

if (NOT TARGET PhysicalModelLib)
    add_subdirectory(../../PhysicalModel PhysicalModelLib)
endif()
if (NOT TARGET CommonLib)
    add_subdirectory(../../Common CommonLib)
endif()

add_library(UtilsLib STATIC 
    CompareCsv.cpp
    WilsonUtils.cpp)

set_target_properties(UtilsLib PROPERTIES PUBLIC_HEADER "CompareCsv.h;WilsonUtils.h")

target_include_directories(UtilsLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(UtilsLib PUBLIC PhysicalModelLib)