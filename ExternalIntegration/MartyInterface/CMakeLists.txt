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

if (NOT TARGET CoreLib)
    add_subdirectory(../../Core CoreLib)
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

add_library(MartyInterfaceLib STATIC 
    CodeGenerator.cpp
    GppCompilerStrategy.cpp
    MakeCompilerStrategy.cpp
    MartyInterface.cpp
    TemplateManager.cpp
    NumericalScriptProcessing/Extractor.cpp
    NumericalScriptProcessing/Interpreter.cpp
    NumericalScriptProcessing/SMParamSetter.cpp
    FileNameManager.cpp
    GeneralModelModifier.cpp
    )

set_target_properties(MartyInterfaceLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_target_properties(MartyInterfaceLib PROPERTIES PUBLIC_HEADER "CinematicExtractor.h;FileNameManager.h;ModelModifier.h;GeneralNumModelModifier.h;CodeGenerator.h;CompilerStrategy.h;GppCompilerStrategy;MakeCompilerStrategy.h;GeneralModelModifier.h;MartyInterface.h;TemplateManager.h;NumericalScriptProcessing/IncludeManager.h;NumericalScriptProcessing/LineProcessor.h;NumericalScriptProcessing/ModelWriter.h;NumericalScriptProcessing/ParamWriter.h;NumericalScriptProcessing/Extractor.h;NumericalScriptProcessing/Interpreter.h;NumericalScriptProcessing/IParamSetter.h;NumericalScriptProcessing/SMParamSetter.h")
# target_include_directories(MartyInterfaceLib PRIVATE BSM)

target_include_directories(MartyInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(MartyInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/NumericalScriptProcessing)
target_link_libraries(MartyInterfaceLib PUBLIC CommonLib DataBaseLib CoreLib)
target_link_libraries(MartyInterfaceLib PUBLIC ConfiguredModule)

add_executable(TEST_MARTY main_marty.cpp)

target_link_libraries(TEST_MARTY MartyInterfaceLib)
