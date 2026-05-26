cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Statistic)
add_compile_options(-O3)

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

# option(ENABLE_SANITIZERS "Enable ASan/UBSan" ON)

# if(ENABLE_SANITIZERS)
#   if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
#     add_compile_options(-g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined)
#     add_link_options(-fsanitize=address,undefined)
#   elseif(MSVC)
#     add_compile_options(/Zi /Od)
#     add_link_options(/INCREMENTAL:NO)
#     # MSVC ASan (pas d'UBSan)
#     add_compile_options(/fsanitize=address)
#   endif()
# endif()

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

# include(FetchContent)

# FetchContent_Declare(minuit_cpp
#   GIT_REPOSITORY https://github.com/doughague/minuit-cpp.git
#   GIT_TAG b5f2558
#   GIT_SHALLOW TRUE
# )

# FetchContent_MakeAvailable(minuit_cpp)

# if(NOT TARGET MinuitCpp)
#   add_library(MinuitCpp INTERFACE)

#   if(TARGET minuit-cpp)
#     target_link_libraries(MinuitCpp INTERFACE minuit-cpp)          # shared
#   elseif(TARGET minuit-cpp-static)
#     target_link_libraries(MinuitCpp INTERFACE minuit-cpp-static)   # static
#   else()
#     message(FATAL_ERROR "Targets minuit-cpp introuvables (minuit-cpp / minuit-cpp-static).")
#   endif()

#   target_include_directories(MinuitCpp INTERFACE
#     ${minuit_cpp_SOURCE_DIR}/include
#   )

#   add_library(MinuitCpp::MinuitCpp ALIAS MinuitCpp)
# endif()

add_library(StatisticLib STATIC 
    adapters/StatCorrelationProxy.cpp
    adapters/StatParameterProxy.cpp
    adapters/StatParamSourcesProxy.cpp
    adapters/StatParamOptimizerProxy.cpp
    adapters/StatDependencyPruner.cpp
    adapters/StatisticInterface.cpp
    adapters/ObservableInterfaceProxy.cpp
    domain/NuisanceReader.cpp
    # domain/RNGHelper.cpp
    domain/JointDistribution.cpp
    domain/Fit.cpp
    # domain/Likelihood.cpp
    domain/CovarianceTransformer.cpp
    domain/Marginal/MarginalConfigFactory.cpp
    domain/MCEngine.cpp
    domain/RvgNuisanceSampler.cpp
    domain/StatisticManager.cpp
    domain/Marginal/MarginalFactory.cpp
    domain/Marginal/GaussianMarginal.cpp
    domain/Marginal/SplitGaussianMarginal.cpp
    domain/Marginal/LikelihoodMarginal.cpp
    domain/Marginal/FlatMarginal.cpp
    domain/Copula/CopulaFactory.cpp
    domain/Copula/GaussianCopula.cpp
    domain/Copula/GenericCopula.cpp
    domain/Copula/StudentTCopula.cpp
    domain/Likelihood/BaseLikelihood.cpp
    domain/Likelihood/IProfilingStrategy.cpp
    domain/Likelihood/ChiSquaredLikelihood.cpp
    domain/Likelihood/ProfiledLikelihood2D.cpp
    domain/Likelihood/Profiler.cpp
    domain/Likelihood/WithGaussianConstraints.cpp
    domain/Contour/AMSContourExtractor.cpp
    domain/Contour/ContourEngine.cpp
    domain/Contour/MnContourExtractor.cpp
)
set_target_properties(StatisticLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(StatisticLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Marginal>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Copula>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Likelihood>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Contour>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain/Marginal>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain/Copula>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain/Likelihood>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Statistic/domain/Contour>
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
target_compile_options(MANAGER_STAT PUBLIC -O3)

add_executable(COPULA_TEST app/test_copulas.cpp)
target_link_libraries(COPULA_TEST StatisticLib)

add_executable(MINUIT_TEST app/test_minuit.cpp)
target_link_libraries(MINUIT_TEST StatisticLib)

add_executable(MINUIT_PIPELINE_TEST app/test_minuit_pipeline.cpp)
target_link_libraries(MINUIT_PIPELINE_TEST StatisticLib)

add_executable(CONTOUR_TEST app/test_contour_minuit.cpp)
target_link_libraries(CONTOUR_TEST StatisticLib)

add_executable(CONTOUR_TEST2 app/test_contour_minuit2.cpp)
target_link_libraries(CONTOUR_TEST2 StatisticLib)

add_executable(REFACTORED_MINUIT app/test_contour_minuit_refactored.cpp)
target_link_libraries(REFACTORED_MINUIT StatisticLib)

add_executable(GAUSSIAN_CONTOUR app/test_gaussian_contour.cpp)
target_link_libraries(GAUSSIAN_CONTOUR StatisticLib)

add_executable(PROFILING_TEST app/test_profiling.cpp)
target_link_libraries(PROFILING_TEST StatisticLib)

add_executable(JPP app/jpp.cpp)
target_link_libraries(JPP StatisticLib)

add_executable(READER app/test_reader.cpp)
target_link_libraries(READER StatisticLib)

add_executable(BKSTARMUMU app/test_contour_BKstarmumu.cpp)
target_link_libraries(BKSTARMUMU StatisticLib)

add_executable(FL app/test_uncertainties_FL.cpp)
target_link_libraries(FL StatisticLib)

add_executable(LIKELIHOODEXP app/likelihood_exploration.cpp)
target_link_libraries(LIKELIHOODEXP StatisticLib)

add_executable(NICELIKELIHOOD app/test_contour_method.cpp)
target_link_libraries(NICELIKELIHOOD StatisticLib)

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