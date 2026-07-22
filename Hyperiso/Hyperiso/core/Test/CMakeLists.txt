cmake_minimum_required(VERSION 3.20)
project(TestHyperiso LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
option(HYPERISO_BUILD_UNIT_TESTS "Build unit tests" ON)
option(HYPERISO_BUILD_INTEGRATION_TESTS "Build integration tests" ON)
option(HYPERISO_BUILD_MARTY_TESTS "Build MartyInterface tests" ON)
option(HYPERISO_BUILD_LEGACY_TESTS "Build legacy/manual smoke tests kept for compatibility" OFF)
option(HYPERISO_BUILD_COMPARISON_TESTS "Build heavy SuperIso/Wilson/observable comparison tests" OFF)

if(NOT TARGET StatisticLib)
    add_subdirectory(../src/Statistic StatisticLib)
endif()
if(NOT TARGET UtilsLib)
    add_subdirectory(Utils UtilsLib)
endif()

include(CTest)
enable_testing()

set(HYPERISO_TEST_ENV "DISABLE_LOGGER=1")

function(hyperiso_add_test target label)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs SOURCES LIBS INCLUDES)
    cmake_parse_arguments(HAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT HAT_SOURCES)
        message(FATAL_ERROR "hyperiso_add_test(${target}) requires SOURCES")
    endif()

    add_executable(${target} ${HAT_SOURCES})
    if(HAT_LIBS)
        target_link_libraries(${target} PRIVATE ${HAT_LIBS})
    endif()
    if(HAT_INCLUDES)
        target_include_directories(${target} PRIVATE ${HAT_INCLUDES})
    endif()

    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES
        ENVIRONMENT "${HYPERISO_TEST_ENV}"
        LABELS "${label}"
    )
endfunction()

function(hyperiso_add_unit target label source)
    if(HYPERISO_BUILD_UNIT_TESTS)
        hyperiso_add_test(${target} "${label};unit" SOURCES ${source} ${ARGN})
    endif()
endfunction()

function(hyperiso_add_integration target label source)
    if(HYPERISO_BUILD_INTEGRATION_TESTS)
        hyperiso_add_test(${target} "${label};integration" SOURCES ${source} ${ARGN})
    endif()
endfunction()

# -----------------------------------------------------------------------------
# Math
# -----------------------------------------------------------------------------
hyperiso_add_unit(testScalarUnit math Math/unit/test_scalar_unit.cpp LIBS MathLib UtilsLib)
hyperiso_add_integration(testScalarIntegration math Math/integration/test_scalar_integration.cpp LIBS MathLib UtilsLib)
hyperiso_add_unit(testMatrixUnit math Math/unit/test_matrix_unit.cpp LIBS MathLib UtilsLib)
hyperiso_add_integration(testMatrixIntegration math Math/integration/test_matrix_integration.cpp LIBS MathLib UtilsLib)

# -----------------------------------------------------------------------------
# Common
# -----------------------------------------------------------------------------
hyperiso_add_unit(testCommon common Common/testLogger.cpp LIBS CommonLib)
hyperiso_add_unit(testlhaidunit common Common/unit/test_lhaid_unit.cpp LIBS CommonLib UtilsLib)
hyperiso_add_integration(testlhaidintegration common Common/integration/test_lhaid_integration.cpp LIBS CommonLib UtilsLib)
hyperiso_add_unit(testblocknameunit common Common/unit/test_blockname_unit.cpp LIBS CommonLib UtilsLib)
hyperiso_add_integration(testblocknameintegration common Common/integration/test_blockname_integration.cpp LIBS CommonLib UtilsLib)
hyperiso_add_unit(testparamidunit common Common/unit/test_paramid_unit.cpp LIBS CommonLib UtilsLib)
hyperiso_add_integration(testparamidintegration common Common/integration/test_paramid_integration.cpp LIBS CommonLib UtilsLib)
hyperiso_add_unit(testenummapperunit common Common/unit/test_enummapper_unit.cpp LIBS CommonLib UtilsLib)
hyperiso_add_integration(testenummapperintegration common Common/integration/test_enummapper_integration.cpp LIBS CommonLib UtilsLib)
hyperiso_add_unit(testDynamicMappersUnit common Common/unit/test_dynamic_mappers_unit.cpp LIBS CommonLib UtilsLib)
hyperiso_add_unit(testObservableTypeIdUnit common Common/unit/test_observable_type_id_unit.cpp LIBS CommonLib UtilsLib)

# -----------------------------------------------------------------------------
# DataBase and Pandas
# -----------------------------------------------------------------------------
hyperiso_add_unit(testlhaelementunit database DataBase/unit/test_lha_elements_unit.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_integration(testlhaelementintegration database DataBase/integration/test_lha_elements_integration.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_unit(testdbnodeunit database DataBase/unit/test_DBnode.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_integration(testdbnodeintegration database DataBase/integration/test_DBnode.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_unit(testdbmanagerunit database DataBase/unit/test_DBmanager.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_integration(testdbmanagerintegration database DataBase/integration/test_DBmanager.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_unit(testlhablocksunit database DataBase/unit/test_lha_blocks_unit.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_integration(testlhablocksintegration database DataBase/integration/test_lha_blocks_integration.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_unit(testjsonparserunit database DataBase/unit/test_jsonparser_unit.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_integration(testjsonparserintegration database DataBase/integration/test_jsonparser_integration.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_unit(testlhaparserunit database DataBase/unit/test_lhaparser_unit.cpp LIBS DataBaseLib UtilsLib)
hyperiso_add_integration(testlhaparserintegration database DataBase/integration/test_lhaparser_integration.cpp LIBS DataBaseLib UtilsLib)

hyperiso_add_unit(testSerieUnit pandas DataBase/Pandas/unit/test_serie_unit.cpp LIBS PandasLib)
hyperiso_add_integration(testSerieIntegration pandas DataBase/Pandas/integration/test_serie_integration.cpp LIBS PandasLib)
hyperiso_add_unit(testDataFrameUnit pandas DataBase/Pandas/unit/test_dataframe_unit.cpp LIBS PandasLib)
hyperiso_add_integration(testDataFrameIntegration pandas DataBase/Pandas/integration/test_dataframe_integration.cpp LIBS PandasLib)

# -----------------------------------------------------------------------------
# Core
# -----------------------------------------------------------------------------
hyperiso_add_unit(testParameterUnit core Core/unit/test_parameter_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_integration(testParameterIntegration core Core/integration/test_parameter_integration.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testBlockUnit core Core/unit/test_block_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_integration(testBlockIntegration core Core/integration/test_block_integration.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testBlockAccessorUnit core Core/unit/test_blockaccessor_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_integration(testBlockAccessorIntegration core Core/integration/test_blockaccessor_integration.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testCoreFileWriterUnit core Core/unit/test_core_file_writer_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testCorrelationRepoUnit core Core/unit/test_correlationrepo_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_integration(testCorrelationRepoIntegration core Core/integration/test_correlationrepo_integration.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testParameterRouterUnit core Core/unit/test_parameterrouter_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_integration(testParameterRouterIntegration core Core/integration/test_parameterrouter_integration.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testDBMementoUnit core Core/unit/test_dbmemento_unit.cpp LIBS CoreLib UtilsLib)
hyperiso_add_integration(testDBMementoIntegration core Core/integration/test_dbmemento_integration.cpp LIBS CoreLib UtilsLib)
hyperiso_add_unit(testMemoryManagerUnit core Core/unit/test_memory_manager_unit.cpp LIBS CoreLib UtilsLib INCLUDES Core/support)
hyperiso_add_integration(testMemoryManagerIntegration core Core/integration/test_memory_manager_integration.cpp LIBS CoreLib UtilsLib INCLUDES Core/support)
hyperiso_add_unit(testParametersUnit core Core/unit/test_parameters_unit.cpp LIBS CoreLib UtilsLib INCLUDES Core/support)
hyperiso_add_integration(testParametersIntegration core Core/integration/test_parameters_integration.cpp LIBS CoreLib UtilsLib INCLUDES Core/support)

# -----------------------------------------------------------------------------
# ExternalIntegration / MartyInterface
# -----------------------------------------------------------------------------
if(HYPERISO_BUILD_MARTY_TESTS)
    hyperiso_add_unit(testCompilerStrategyUnit martyinterface ExternalIntegration/MartyInterface/unit/test_compiler_strategy_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_integration(testCompilerStrategyIntegration martyinterface ExternalIntegration/MartyInterface/integration/test_compiler_strategy_integration.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testModelFileChecker martyinterface ExternalIntegration/MartyInterface/unit/test_model_file_checker_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testGeneralModelModifierUnit martyinterface ExternalIntegration/MartyInterface/unit/test_general_model_modifier_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testTemplateManagerUnit martyinterface ExternalIntegration/MartyInterface/unit/test_template_manager_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_integration(testGeneralModelModifierIntegration martyinterface ExternalIntegration/MartyInterface/integration/test_general_model_modifier_integration.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testExtractorUnit martyinterface ExternalIntegration/MartyInterface/unit/test_extractor_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testCinematicExtractorUnit martyinterface ExternalIntegration/MartyInterface/unit/test_cinematic_extractor_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testFileWriterUnit martyinterface ExternalIntegration/MartyInterface/unit/test_file_writer_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testLineProcessorUnit martyinterface ExternalIntegration/MartyInterface/unit/test_line_processor_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testParamWriterUnit martyinterface ExternalIntegration/MartyInterface/unit/test_param_writer_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_unit(testSMParamSetterUnit martyinterface ExternalIntegration/MartyInterface/unit/test_param_setter_unit.cpp LIBS MartyInterfaceLib UtilsLib)
    hyperiso_add_integration(testWriterPipelineIntegration martyinterface ExternalIntegration/MartyInterface/integration/test_writer_pipeline_integration.cpp LIBS MartyInterfaceLib UtilsLib)
endif()

# -----------------------------------------------------------------------------
# PhysicalModel
# -----------------------------------------------------------------------------
hyperiso_add_unit(testWilsonUnit physicalmodel PhysicalModel/unit/test_wilson_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testWilsonIntegration physicalmodel PhysicalModel/integration/test_wilson_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testMartyWilsonUnit physicalmodel PhysicalModel/unit/test_marty_wilson_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testMartyWilsonIntegration physicalmodel PhysicalModel/integration/test_marty_wilson_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testWilsonGroupUnit physicalmodel PhysicalModel/unit/test_wilsongroup_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testWilsonGroupIntegration physicalmodel PhysicalModel/integration/test_wilsongroup_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testWilsonRegistryUnit physicalmodel PhysicalModel/unit/test_wilson_registry_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testWilsonRegistryIntegration physicalmodel PhysicalModel/integration/test_wilson_registry_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testCoefficientGroupBuilderUnit physicalmodel PhysicalModel/unit/test_coefficient_group_builder_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testCoefficientGroupBuilderIntegration physicalmodel PhysicalModel/integration/test_coefficient_group_builder_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testGroupDefinitionUnit physicalmodel PhysicalModel/unit/test_group_definition_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testGroupDefinitionIntegration physicalmodel PhysicalModel/integration/test_group_definition_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testCustomWilsonUnit physicalmodel PhysicalModel/unit/test_custom_wilson_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testCustomWilsonIntegration physicalmodel PhysicalModel/integration/test_custom_wilson_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testCustomWilsonLambdaConfigUnit physicalmodel PhysicalModel/unit/test_custom_wilson_lambda_config_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testWilsonManagerUnit physicalmodel PhysicalModel/unit/test_wilson_manager_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testWilsonManagerIntegration physicalmodel PhysicalModel/integration/test_wilson_manager_integration.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_unit(testHasWilsonUnit physicalmodel PhysicalModel/unit/test_has_wilson_unit.cpp LIBS PhysicalModelLib UtilsLib)
hyperiso_add_integration(testHasWilsonIntegration physicalmodel PhysicalModel/integration/test_has_wilson_integration.cpp LIBS PhysicalModelLib UtilsLib)

# -----------------------------------------------------------------------------
# BusinessLogic
# -----------------------------------------------------------------------------
hyperiso_add_unit(testObsWilsonHelperUnit businesslogic BusinessLogic/unit/test_obswilsonhelper_unit.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_integration(testObsWilsonHelperIntegration businesslogic BusinessLogic/integration/test_obswilsonhelper_integration.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_unit(testDecayParentUnit businesslogic BusinessLogic/unit/test_decayparent_unit.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_integration(testDecayParentIntegration businesslogic BusinessLogic/integration/test_decayparent_integration.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_unit(testObservableUnit businesslogic BusinessLogic/unit/test_observable_unit.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_integration(testObservableIntegration businesslogic BusinessLogic/integration/test_observable_integration.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_unit(testObsManagerUnit businesslogic BusinessLogic/unit/test_obsmanager_unit.cpp LIBS BusinessLogicLib UtilsLib)
hyperiso_add_integration(testObsManagerIntegration businesslogic BusinessLogic/integration/test_obsmanager_integration.cpp LIBS BusinessLogicLib UtilsLib)

# -----------------------------------------------------------------------------
# Statistic
# -----------------------------------------------------------------------------
hyperiso_add_unit(testRuntimeNuisancePathsAdapterUnit statistic Statistic/unit/test_runtime_nuisance_paths_adapter_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testMarginalsUnit statistic Statistic/unit/test_marginals_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testMarginalsIntegration statistic Statistic/integration/test_marginals_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testCopulasUnit statistic Statistic/unit/test_copulas_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testCopulasIntegration statistic Statistic/integration/test_copulas_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testJointDistributionUnit statistic Statistic/unit/test_jointdistribution_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testJointDistributionIntegration statistic Statistic/integration/test_jointdistribution_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testLikelihoodBaseUnit statistic Statistic/unit/test_likelihood_base_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testLikelihoodBaseIntegration statistic Statistic/integration/test_likelihood_base_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testGradientHelperUnit statistic Statistic/unit/test_gradienthelper_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testGradientHelperIntegration statistic Statistic/integration/test_gradienthelper_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testProfilerUnit statistic Statistic/unit/test_profiler_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testProfilerIntegration statistic Statistic/integration/test_profiler_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testIndexingUnit statistic Statistic/unit/test_indexing_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testIndexingIntegration statistic Statistic/integration/test_indexing_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testChiSquaredLikelihoodUnit statistic Statistic/unit/test_chisquaredlikelihood_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testChiSquaredLikelihoodIntegration statistic Statistic/integration/test_chisquaredlikelihood_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testFitProfileHessianUnit statistic Statistic/unit/test_fit_profile_hessian_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testProfilingStrategyUnit statistic Statistic/unit/test_profilingstrategy_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testProfilingStrategyIntegration statistic Statistic/integration/test_profilingstrategy_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testProfiledLikelihood2DUnit statistic Statistic/unit/test_profiledlikelihood2d_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testProfiledLikelihood2DIntegration statistic Statistic/integration/test_profiledlikelihood2d_integration.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_unit(testWithGaussianConstraintsUnit statistic Statistic/unit/test_withgaussianconstraints_unit.cpp LIBS StatisticLib UtilsLib)
hyperiso_add_integration(testWithGaussianConstraintsIntegration statistic Statistic/integration/test_withgaussianconstraints_integration.cpp LIBS StatisticLib UtilsLib)

# -----------------------------------------------------------------------------
# Legacy/manual smoke tests. Kept out of the default build because several of
# them depend on old comparison files, absolute environments, or optional BSM
# backends. Enable explicitly when maintaining those workflows.
# -----------------------------------------------------------------------------
if(HYPERISO_BUILD_LEGACY_TESTS)
    hyperiso_add_test(testReader "database;legacy" SOURCES DataBase/testReader.cpp LIBS DataBaseLib ConfiguredModule)
    hyperiso_add_test(testjsonreader "database;legacy" SOURCES DataBase/testJsonReader.cpp LIBS DataBaseLib ConfiguredModule)
    hyperiso_add_test(testmappingjson "database;legacy" SOURCES DataBase/testMapping.cpp LIBS DataBaseLib ConfiguredModule)
    hyperiso_add_test(teststatjson "database;legacy" SOURCES DataBase/testJsonStatParser.cpp LIBS DataBaseLib ConfiguredModule)

    hyperiso_add_test(testSMparams "core;legacy" SOURCES Core/testSMParameters.cpp LIBS CoreLib)
    hyperiso_add_test(testFlavorParams "core;legacy" SOURCES Core/testFlavorParams.cpp LIBS CoreLib)
    hyperiso_add_test(testObservableLegacy "businesslogic;legacy" SOURCES BusinessLogic/testObservable.cpp LIBS BusinessLogicLib)
    hyperiso_add_test(testChi2Legacy "businesslogic;legacy" SOURCES BusinessLogic/testChi2.cpp LIBS BusinessLogicLib)
    hyperiso_add_test(testCustomWilsonLegacy "businesslogic;legacy" SOURCES BusinessLogic/testCustomWilson.cpp LIBS BusinessLogicLib)
endif()

# -----------------------------------------------------------------------------
# Heavy numerical comparison tests.
# -----------------------------------------------------------------------------
if(HYPERISO_BUILD_COMPARISON_TESTS)
    hyperiso_add_test(testSMWilson "wilson;comparison" SOURCES Wilson/testSMWilson.cpp LIBS PhysicalModelLib UtilsLib)
    hyperiso_add_test(testSMCQPrimeWilson "wilson;comparison" SOURCES Wilson/testSMCQPrimeWilson.cpp LIBS PhysicalModelLib UtilsLib)
    hyperiso_add_test(testBase1Wilson "wilson;comparison" SOURCES Wilson/testBase1Wilson.cpp LIBS PhysicalModelLib UtilsLib)
    hyperiso_add_test(testBase2Wilson "wilson;comparison" SOURCES Wilson/testBase2Wilson.cpp LIBS PhysicalModelLib UtilsLib)
    # hyperiso_add_test(testObservables "observables;comparison" SOURCES Observables/test_observables.cpp LIBS BusinessLogicLib UtilsLib)

    hyperiso_add_test(testTHDMWilson "wilson;comparison;2hdmc" SOURCES Wilson/testTHDMWilson.cpp LIBS PhysicalModelLib UtilsLib)
    hyperiso_add_test(testTHDMCQPrimeWilson "wilson;comparison;2hdmc" SOURCES Wilson/testTHDMCQPrimeWilson.cpp LIBS PhysicalModelLib UtilsLib)

    if(BUILD_WITH_SOFTSUSY)
        hyperiso_add_test(testSUSYWilson "wilson;comparison;softsusy" SOURCES Wilson/testSUSYWilson.cpp LIBS PhysicalModelLib UtilsLib)
        hyperiso_add_test(testSUSYCQPrimeWilson "wilson;comparison;softsusy" SOURCES Wilson/testSUSYCQPrimeWilson.cpp LIBS PhysicalModelLib UtilsLib)
    endif()
endif()
