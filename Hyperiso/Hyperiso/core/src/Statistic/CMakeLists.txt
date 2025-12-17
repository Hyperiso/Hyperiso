cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Statistic)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
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

add_library(StatisticLib STATIC 
    adapters/StatCorrelationProxy.cpp
    adapters/StatParameterProxy.cpp
    adapters/StatParamSourcesProxy.cpp
    adapters/StatParamOptimizerProxy.cpp
    domain/RNGHelper.cpp
    domain/RandomVectorGenerator.cpp
    domain/Fit.cpp
    domain/Likelihood.cpp
    domain/CovarianceTransformer.cpp
    domain/Distribution/DistributionFactory.cpp
    domain/Distribution/StandardNormal.cpp
    domain/Distribution/LikelihoodDiscrete.cpp
    domain/Distribution/StandardFlat.cpp
)
set_target_properties(StatisticLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(StatisticLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Distribution>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain/Distribution>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/ports>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/adapters>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(StatisticLib PUBLIC BusinessLogicLib PhysicalModelLib CoreLib CommonLib ConfiguredModule MathLib)

add_executable(RNG_TEST app/main_rng.cpp)
target_link_libraries(RNG_TEST StatisticLib)

add_executable(STAT_TEST app/example_main.cpp)
target_link_libraries(STAT_TEST StatisticLib)

add_executable(REAL_STAT_TEST app/example_main_adapter.cpp)
target_link_libraries(REAL_STAT_TEST StatisticLib)

add_executable(ADAPTERS_TEST app/test_adapters.cpp)
target_link_libraries(ADAPTERS_TEST StatisticLib)

add_executable(MANAGER_STAT app/example_manager.cpp)
target_link_libraries(MANAGER_STAT StatisticLib)

include(GNUInstallDirs)

install(TARGETS StatisticLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/StatisticLib
  FILES_MATCHING
    PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)