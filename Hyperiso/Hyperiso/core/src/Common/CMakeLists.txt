cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Common)

add_library(CommonLib STATIC 
    Logger.cpp
    General.cpp
    Observer.cpp
    Mapper/Map.cpp)

set_target_properties(CommonLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(CommonLib PROPERTIES PUBLIC_HEADER "Logger.h;General.h;Utils.h")

target_include_directories(CommonLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(CommonLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Mapper)


add_executable(aaaah aaaah.cpp)

target_link_libraries(aaaah CommonLib)