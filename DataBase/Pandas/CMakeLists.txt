cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Pandas)

if (NOT TARGET CommonLib)
    add_subdirectory(../../Common CommonLib)
endif()

add_library(PandasLib STATIC
    DataFrame.cpp
    CSVReader.cpp
    )

set_target_properties(PandasLib PROPERTIES PUBLIC_HEADER "Series.h;DataFrame.h;CSVReader.h;CSVOptions.h")

target_include_directories(PandasLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(PandasLib PUBLIC CommonLib)

add_executable(TEST_PANDAS main.cpp)

target_link_libraries(TEST_PANDAS PRIVATE PandasLib)