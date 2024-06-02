cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Core)

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()
if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured ConfiguredModule)
endif()
if (NOT TARGET DataBaseLib)
    add_subdirectory(../DataBase DataBaseLib)
endif()
if (NOT TARGET MathLib)
    add_subdirectory(../Math MathLib)
endif()
if (NOT TARGET ExternalIntegrationLib)
    add_subdirectory(../ExternalIntegration ExternalIntegrationLib)
endif()


add_library(CoreLib STATIC 
    MemoryManager.cpp
    Parameters.cpp
    QCDParameters.cpp)

set_target_properties(CoreLib PROPERTIES PUBLIC_HEADER "MemoryManager.h;Parameters.h;QCDParameters.h")

target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(CoreLib PUBLIC CommonLib DataBaseLib ConfiguredLib ExternalIntegrationLib MathLib)