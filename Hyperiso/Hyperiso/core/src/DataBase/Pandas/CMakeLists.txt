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

set_target_properties(PandasLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
# set_target_properties(PandasLib PROPERTIES PUBLIC_HEADER "Series.h;DataFrame.h;CSVReader.h;CSVOptions.h")

# target_include_directories(PandasLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


target_include_directories(PandasLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase/Pandas>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(PandasLib PUBLIC CommonLib)

include(GNUInstallDirs)

install(TARGETS PandasLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/DataBase/Pandas
  FILES_MATCHING
    PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)

add_executable(TEST_PANDAS main.cpp)

target_link_libraries(TEST_PANDAS PRIVATE PandasLib)