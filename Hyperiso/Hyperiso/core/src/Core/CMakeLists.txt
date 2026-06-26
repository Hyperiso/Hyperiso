cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)



project(Core)

if (NOT DEFINED BUILD_WITH_APP)
    option(BUILD_WITH_APP "Build app/ executable targets" OFF)
endif()

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MINUIT2)
    option(BUILD_WITH_MINUIT2 "Build project with Minuit2 support" OFF)
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
    adapters/DBMemento.cpp
    adapters/ParamBlockLoader.cpp
    adapters/ParamBlockWriter.cpp
    adapters/FileWriter.cpp
    adapters/APIAdapter.cpp
    adapters/MartyAdapter.cpp
    adapters/HyperisoMaster.cpp
    adapters/CorrelationProvider.cpp
    adapters/ParameterProvider.cpp
    adapters/SpectrumCalculator.cpp
    adapters/CorrelationAdapter.cpp
    adapters/CompositeParamCreator.cpp
    adapters/DependantBlockInfoProvider.cpp
    adapters/QCDProvider.cpp
    adapters/QEDProvider.cpp
    adapters/ParameterSetter.cpp
    adapters/ParameterShifter.cpp
    adapters/Freezer.cpp
    adapters/BlockProvider.cpp
    adapters/ParamSourcesProvider.cpp
    adapters/ParamOptimizerAdapter.cpp
    adapters/DependencyPruner.cpp
    domain/CorrelationRepo.cpp
    domain/ParameterRouter.cpp
    domain/Parameters.cpp
    domain/Parameter.cpp
    domain/DependentParameter.cpp
    domain/MemoryManager.cpp
    domain/ParameterRuntimeContext.cpp
    domain/Block.cpp
    domain/BlockAccessor.cpp
    domain/SourceView.cpp
    domain/DependentBlockManager.cpp
    domain/ParamOptimizer.cpp
    infrastructure/QCDHelper.cpp
    infrastructure/EWHelper.cpp
    )
set_target_properties(CoreLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# set_target_properties(CoreLib PROPERTIES PUBLIC_HEADER "MemoryManager.h;Parameters.h;Parameter.h;QCDHelper.h;BlocksCreator.h;BlockAccessor.h;CorrelationRepo.h;CorrelationCreator.h;ParameterRouter.h;")

# target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
# target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/domain)
# target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/infrastructure)
# target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ports)



target_include_directories(CoreLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/infrastructure>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>

    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Core>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Core/adapters>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Core/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Core/infrastructure>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Core/ports>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(CoreLib PUBLIC CommonLib DataBaseLib ConfiguredModule ExternalIntegrationLib MathLib)

if (${ENABLE_WARNINGS})
    target_set_warnings(CoreLib
                        ${ENABLE_WARNINGS}
                        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if (BUILD_WITH_APP)
add_executable(TEST_CORE app/main_core.cpp)
target_link_libraries(TEST_CORE CoreLib)

add_executable(OPTIMIZER_TEST app/test_optimizer.cpp)
target_link_libraries(OPTIMIZER_TEST CoreLib)

add_executable(SMALL_OPTIMIZER_TEST app/test_little_optimizer.cpp)
target_link_libraries(SMALL_OPTIMIZER_TEST CoreLib)

add_executable(DEEP_OPTIMIZER_TEST app/main_optimizer_deeper.cpp)
target_link_libraries(DEEP_OPTIMIZER_TEST CoreLib)

add_executable(PRINT_OPTIMIZER_TEST app/main_optimizer_print.cpp)
target_link_libraries(PRINT_OPTIMIZER_TEST CoreLib)

add_executable(TEST_MMH app/main_parameters.cpp)
target_link_libraries(TEST_MMH CoreLib)

add_executable(LAZY_TEST app/main_lazy.cpp)
target_link_libraries(LAZY_TEST CoreLib)

add_executable(LAZY_TEST_2 app/main_more_lazy.cpp)
target_link_libraries(LAZY_TEST_2 CoreLib)

add_executable(MORE_SPEED app/main_speed.cpp)
target_link_libraries(MORE_SPEED CoreLib)

add_executable(DETACH app/main_detach.cpp)
target_link_libraries(DETACH CoreLib)
endif()

include(GNUInstallDirs)

install(TARGETS CoreLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Core
        FILES_MATCHING
          PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
          PATTERN "*.cpp" EXCLUDE
          PATTERN "app" EXCLUDE
          PATTERN "CMakeFiles" EXCLUDE)