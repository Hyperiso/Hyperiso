cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(ExternalMartyIntegration)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()

if (NOT TARGET CommonLib)
    add_subdirectory(../../Common CommonLib)
endif()

if (NOT TARGET ExternalIntegrationLib)
    add_subdirectory(.. ExternalIntegrationLib)
endif()

if (NOT TARGET DataBaseLib)
    add_subdirectory(../../DataBase DataBaseLib)
endif()

if (NOT TARGET ConfiguredModule)
    add_subdirectory(../../Configured Configured)
endif()

add_library(ExternalMartyIntegrationLib STATIC 
    CodeGenerator.cpp
    GppCompilerStrategy.cpp
    MakeCompilerStrategy.cpp
    MartyInterface.cpp
    SMModelModifier.cpp
    TemplateManager.cpp
    Extractor.cpp
    Interpreter.cpp
    SMParamSetter.cpp
    )

set_target_properties(ExternalMartyIntegrationLib PROPERTIES PUBLIC_HEADER "ModelModifier.h;GeneralNumModelModifier.h;NumModelModifier.h;CodeGenerator.h;CompilerStrategy.h;GppCompilerStrategy;MakeCompilerStrategy.h;MartyInterface.h;SMModelModifier.h;TemplateManager.h;IncludeManager.h;LineProcessor.h;ModelWriter.h;ParamWriter.h;Extractor.h;Interpreter.h;IParamSetter.h;SMParamSetter.h")
# target_include_directories(ExternalMartyIntegrationLib PRIVATE BSM)
target_link_libraries(ExternalMartyIntegrationLib PUBLIC CommonLib DataBaseLib)
target_link_libraries(ExternalMartyIntegrationLib PUBLIC ConfiguredModule)

add_executable(TEST_MARTY main_marty.cpp)

target_link_libraries(TEST_MARTY ExternalMartyIntegrationLib)
