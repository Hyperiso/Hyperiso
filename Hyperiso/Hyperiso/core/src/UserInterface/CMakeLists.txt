cmake_minimum_required(VERSION 3.20)

project(UserInterface LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MINUIT2)
    option(BUILD_WITH_MINUIT2 "Build project with Minuit2 support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_APP)
    option(BUILD_WITH_APP "Build Hyperiso command-line applications" OFF)
endif()

if (NOT TARGET StatisticLib)
    add_subdirectory(../Statistic StatisticLib)
endif()
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

add_library(UserInterfaceLib STATIC
    CliUtils.cpp
    WilsonHandler.cpp
    ObservableHandler.cpp
    StatisticHandler.cpp
    ArgsParser.cpp
    YamlInputReader.cpp
    UserParameterProxy.cpp
    CsvWriter.cpp
    JsonWriter.cpp
    ScanRunner.cpp
    WilsonExtractor.cpp
    ObservableExtractor.cpp
)

set_target_properties(UserInterfaceLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(UserInterfaceLib
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/UserInterface>
)

target_link_libraries(UserInterfaceLib
  PUBLIC
    StatisticLib
    BusinessLogicLib
    PhysicalModelLib
    CoreLib
    CommonLib
    ConfiguredModule
    MathLib
)

if(DEFINED ENABLE_CLANG_TIDY AND ENABLE_CLANG_TIDY)
    add_clang_tidy_to_target(UserInterfaceLib)
endif()

if(BUILD_WITH_APP)
    add_executable(hyperiso-ui main.cpp)
    target_link_libraries(hyperiso-ui PRIVATE UserInterfaceLib)
endif()

include(GNUInstallDirs)

install(TARGETS UserInterfaceLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(BUILD_WITH_APP)
    install(TARGETS hyperiso-ui
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/UserInterface
  FILES_MATCHING
    PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)
