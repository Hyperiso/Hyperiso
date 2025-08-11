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
    GeneralNumModelModifier.cpp
    CompilerStrategy.cpp
    CinematicExtractor.cpp
    adapters/MartyParameterProxy.cpp
    )

set_target_properties(MartyInterfaceLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_target_properties(MartyInterfaceLib PROPERTIES PUBLIC_HEADER "CinematicExtractor.h;FileNameManager.h;ModelModifier.h;GeneralNumModelModifier.h;CodeGenerator.h;CompilerStrategy.h;GppCompilerStrategy;MakeCompilerStrategy.h;GeneralModelModifier.h;MartyInterface.h;TemplateManager.h;NumericalScriptProcessing/IncludeManager.h;NumericalScriptProcessing/LineProcessor.h;NumericalScriptProcessing/ModelWriter.h;NumericalScriptProcessing/ParamWriter.h;NumericalScriptProcessing/Extractor.h;NumericalScriptProcessing/Interpreter.h;NumericalScriptProcessing/IParamSetter.h;NumericalScriptProcessing/SMParamSetter.h")
# target_include_directories(MartyInterfaceLib PRIVATE BSM)

# target_include_directories(MartyInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_include_directories(MartyInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ports)
# target_include_directories(MartyInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
# target_include_directories(MartyInterfaceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/NumericalScriptProcessing)

# target_include_directories(MartyInterfaceLib
#   PUBLIC
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
#     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
#   PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}
#     ${CMAKE_CURRENT_SOURCE_DIR}/ports
#     ${CMAKE_CURRENT_SOURCE_DIR}/adapters
#     ${CMAKE_CURRENT_SOURCE_DIR}/NumericalScriptProcessing
# )
target_include_directories(MartyInterfaceLib
  PUBLIC
    # build: expose exactement où vivent tes headers inclus sans préfixe
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/NumericalScriptProcessing>

    # install: miroirs sous include/
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration/MartyInterface>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration/MartyInterface/ports>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration/MartyInterface/adapters>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration/MartyInterface/NumericalScriptProcessing>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MartyInterfaceLib PUBLIC CommonLib DataBaseLib CoreLib)
target_link_libraries(MartyInterfaceLib PUBLIC ConfiguredModule)

include(GNUInstallDirs)

install(TARGETS MartyInterfaceLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration/MartyInterface
  FILES_MATCHING
    PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)

add_executable(TEST_MARTY main_marty.cpp)

target_link_libraries(TEST_MARTY MartyInterfaceLib)
