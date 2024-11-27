cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(UserInterface)
add_library(UserInterfaceLib WilsonHandler.cpp ArgsParser.cpp ObservableHandler.cpp)

set_target_properties(UserInterfaceLib PROPERTIES PUBLIC_HEADER "ObservableHandler.h;ArgsParser.h;WilsonHandler.h;HyperisoInterface.h")

target_include_directories(UserInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if (NOT TARGET BusinessLogicLib)
    add_subdirectory(../BusinessLogic BusinessLogicLib)
endif()

if (NOT TARGET PhysicalModelLib)
    add_subdirectory(../PhysicalModel PhysicalModelLib)
endif()

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()
if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured ConfiguredModule)
endif()
if (NOT TARGET MathLib)
    add_subdirectory(../Math MathLib)
endif()
if (NOT TARGET CoreLib)
    add_subdirectory(../Core CoreLib)
endif()

target_link_libraries(UserInterfaceLib PUBLIC BusinessLogicLib PhysicalModelLib CommonLib ConfiguredModule MathLib CoreLib)

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target(UserInterfaceLib)
endif()

add_executable(main_user main.cpp)
target_link_libraries(main_user UserInterfaceLib)