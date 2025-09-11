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
    domain/lha_blocks.cpp
    domain/lha_elements.cpp
    domain/DBNode.cpp
    domain/DBManager.cpp
    Service/MappingDataBase.cpp
    adapters/JsonParser.cpp
    adapters/LhaParser.cpp
    adapters/YamlParser.cpp
    adapters/DictNodeProvider.cpp
    adapters/LhaNodeProvider.cpp
    adapters/NodeProviderFactory.cpp
)

set_target_properties(DataBaseLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Service)
# target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
# target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/domain)
# target_include_directories(DataBaseLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ports)


target_include_directories(DataBaseLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Service>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>
    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase/adapters>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase/Service>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase/ports>
)

target_link_libraries(DataBaseLib PUBLIC PandasLib CommonLib)

# add_executable(TEST_DB app/main_db.cpp)
# target_link_libraries(TEST_DB DataBaseLib)

add_executable(db_test app/main_db.cpp)

target_link_libraries(db_test DataBaseLib)

add_executable(TEST_PARSER app/main_parser.cpp)

target_link_libraries(TEST_PARSER DataBaseLib)

include(GNUInstallDirs)

install(TARGETS DataBaseLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase
  FILES_MATCHING
    PATTERN "*.h"    PATTERN "*.hpp"
    PATTERN "*.tpp"  PATTERN "*.ipp"  PATTERN "*.tcc"  PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "app"  EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)