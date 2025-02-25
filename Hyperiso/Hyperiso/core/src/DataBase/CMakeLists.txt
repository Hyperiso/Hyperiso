cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(DataBase)

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()

if (NOT TARGET PandasLib)
    add_subdirectory(Pandas PandasLib)
endif()

add_library(DataBaseLib STATIC
    Entity/lha_blocks.cpp
    Entity/lha_elements.cpp
    Entity/JsonParameters.cpp
    Entity/DBNode.cpp
    Service/lha_reader.cpp
    Service/MappingDataBase.cpp
    Parser/json_parser.cpp
    Parser/lha_parser.cpp
    Parser/Parser.cpp
)

set_target_properties(DataBaseLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(DataBaseLib PROPERTIES PUBLIC_HEADER "Service/lha_reader.h")

target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Entity)
target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Service)
target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Parser)

target_link_libraries(DataBaseLib PUBLIC PandasLib CommonLib)

add_executable(TEST_DB main_db.cpp)

target_link_libraries(TEST_DB DataBaseLib)

add_executable(TEST_json main_json.cpp)

target_link_libraries(TEST_json DataBaseLib)

add_executable(TEST_PARSER main_parser.cpp)

target_link_libraries(TEST_PARSER DataBaseLib)