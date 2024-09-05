cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Core)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

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

target_link_libraries(CoreLib PUBLIC CommonLib DataBaseLib ConfiguredModule ExternalIntegrationLib MathLib)

if (${ENABLE_WARNINGS})
    target_set_warnings(CoreLib
                        ${ENABLE_WARNINGS}
                        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

# add_executable(TEST main_parameters.cpp)
# target_link_libraries(TEST CoreLib)