cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT TARGET StatisticLib)
    add_subdirectory(../src/Statistic StatisticLib)
endif()
if (NOT TARGET UtilsLib)
    add_subdirectory(Utils UtilsLib)
endif()

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

# include(Catch)
set(CTEST_ENVIRONMENT_VARIABLE "CTEST_RUN")
# add_executable(testSMWilson testSMWilson.cpp)
# target_link_libraries(testSMWilson PhysicalModelLib)

# add_executable(testSusyWilson testSusyWilson.cpp)
# target_link_libraries(testSusyWilson PhysicalModelLib)

# add_executable(testTHDMWilson testTHDMWilson.cpp)
# target_link_libraries(testTHDMWilson PhysicalModelLib)

# add_executable(testWilson testWilson.cpp)
# target_link_libraries(testWilson PhysicalModelLib)

#TODO : Observable tests
# add_executable(testObservable BusinessLogic/testObservable.cpp)
# target_link_libraries(testObservable BusinessLogicLib)

# add_executable(testStat BusinessLogic/testChi2.cpp)
# target_link_libraries(testStat BusinessLogicLib)

# add_executable(testCustomWilsons BusinessLogic/testCustomWilson.cpp)
# target_link_libraries(testCustomWilsons BusinessLogicLib)

# add_executable(testParams testParams.cpp)
# target_link_libraries(testParams CoreLib)

# add_executable(testPrimeCQWilson testPrimeCQWilson.cpp)
# target_link_libraries(testPrimeCQWilson PhysicalModelLib)

# add_executable(testPrimeCQTHDMWilson testPrimeCQTHDMWilson.cpp)
# target_link_libraries(testPrimeCQTHDMWilson PhysicalModelLib)

# add_executable(testPrimeCQSUSYWilson testPrimeCQSUSYWilson.cpp)
# target_link_libraries(testPrimeCQSUSYWilson PhysicalModelLib)



# add_executable(test2HDMC test2HDMC.cpp)
# target_link_libraries(test2HDMC ExternalIntegrationLib)

# add_executable(testQCDParameters testQCDParameters.cpp)
# target_link_libraries(testQCDParameters CoreLib)

# add_executable(testReader testReader.cpp)
# target_link_libraries(testReader DataBaseLib)

add_executable(testScalarIntegration Math/integration/test_scalar_integration.cpp)
target_link_libraries(testScalarIntegration MathLib UtilsLib)

add_executable(testScalarUnit Math/unit/test_scalar_unit.cpp)
target_link_libraries(testScalarUnit MathLib UtilsLib)

add_executable(testMatrixIntegration Math/integration/test_matrix_integration.cpp)
target_link_libraries(testMatrixIntegration MathLib UtilsLib)

add_executable(testMatrixUnit Math/unit/test_matrix_unit.cpp)
target_link_libraries(testMatrixUnit MathLib UtilsLib)

#Common

add_executable(testCommon Common/testLogger.cpp)
target_link_libraries(testCommon CommonLib)

add_executable(testlhaidunit Common/unit/test_lhaid_unit.cpp)
target_link_libraries(testlhaidunit CommonLib UtilsLib)

add_executable(testlhaidintegration Common/integration/test_lhaid_integration.cpp)
target_link_libraries(testlhaidintegration CommonLib UtilsLib)

add_executable(testblocknameunit Common/unit/test_blockname_unit.cpp)
target_link_libraries(testblocknameunit CommonLib UtilsLib)

add_executable(testblocknameintegration Common/integration/test_blockname_integration.cpp)
target_link_libraries(testblocknameintegration CommonLib UtilsLib)

add_executable(testparamidunit Common/unit/test_paramid_unit.cpp)
target_link_libraries(testparamidunit CommonLib UtilsLib)

add_executable(testparamidintegration Common/integration/test_paramid_integration.cpp)
target_link_libraries(testparamidintegration CommonLib UtilsLib)

add_executable(testenummapperunit Common/unit/test_enummapper_unit.cpp)
target_link_libraries(testenummapperunit CommonLib UtilsLib)

add_executable(testenummapperintegration Common/integration/test_enummapper_integration.cpp)
target_link_libraries(testenummapperintegration CommonLib UtilsLib)

#Database

add_executable(testReader DataBase/testReader.cpp)
target_link_libraries(testReader DataBaseLib ConfiguredModule)

add_executable(testjsonreader DataBase/testJsonReader.cpp)
target_link_libraries(testjsonreader DataBaseLib ConfiguredModule)

add_executable(testmappingjson DataBase/testMapping.cpp)
target_link_libraries(testmappingjson DataBaseLib ConfiguredModule)

add_executable(teststatjson DataBase/testJsonStatParser.cpp)
target_link_libraries(teststatjson DataBaseLib ConfiguredModule)

add_executable(testlhaelementunit DataBase/unit/test_lha_elements_unit.cpp)
target_link_libraries(testlhaelementunit DataBaseLib UtilsLib)

add_executable(testlhaelementintegration DataBase/integration/test_lha_elements_integration.cpp)
target_link_libraries(testlhaelementintegration DataBaseLib UtilsLib)

add_executable(testdbnodeunit DataBase/unit/test_DBnode.cpp)
target_link_libraries(testdbnodeunit DataBaseLib UtilsLib)

add_executable(testdbnodeintegration DataBase/integration/test_DBnode.cpp)
target_link_libraries(testdbnodeintegration DataBaseLib UtilsLib)

add_executable(testdbmanagerunit DataBase/unit/test_DBmanager.cpp)
target_link_libraries(testdbmanagerunit DataBaseLib UtilsLib)

add_executable(testdbmanagerintegration DataBase/integration/test_DBmanager.cpp)
target_link_libraries(testdbmanagerintegration DataBaseLib UtilsLib)

add_executable(testlhablocksunit DataBase/unit/test_lha_blocks_unit.cpp)
target_link_libraries(testlhablocksunit DataBaseLib UtilsLib)

add_executable(testlhablocksintegration DataBase/integration/test_lha_blocks_integration.cpp)
target_link_libraries(testlhablocksintegration DataBaseLib UtilsLib)

add_executable(testjsonparserunit DataBase/unit/test_jsonparser_unit.cpp)
target_link_libraries(testjsonparserunit DataBaseLib UtilsLib)

add_executable(testjsonparserintegration DataBase/integration/test_jsonparser_integration.cpp)
target_link_libraries(testjsonparserintegration DataBaseLib UtilsLib)

add_executable(testlhaparserunit DataBase/unit/test_lhaparser_unit.cpp)
target_link_libraries(testlhaparserunit DataBaseLib UtilsLib)

add_executable(testlhaparserintegration DataBase/integration/test_lhaparser_integration.cpp)
target_link_libraries(testlhaparserintegration DataBaseLib UtilsLib)

#Pandas

add_executable(testSerieUnit DataBase/Pandas/unit/test_serie_unit.cpp)
target_link_libraries(testSerieUnit PandasLib)

add_executable(testSerieIntegration DataBase/Pandas/integration/test_serie_integration.cpp)
target_link_libraries(testSerieIntegration PandasLib)

add_executable(testDataFrameUnit DataBase/Pandas/unit/test_dataframe_unit.cpp)
target_link_libraries(testDataFrameUnit PandasLib)

add_executable(testDataFrameIntegration DataBase/Pandas/integration/test_dataframe_integration.cpp)
target_link_libraries(testDataFrameIntegration PandasLib)

#Core


add_executable(testParameterIntegration Core/integration/test_parameter_integration.cpp)
target_link_libraries(testParameterIntegration CoreLib UtilsLib)

add_executable(testParameterUnit Core/unit/test_parameter_unit.cpp)
target_link_libraries(testParameterUnit CoreLib UtilsLib)

add_executable(testBlockIntegration Core/integration/test_block_integration.cpp)
target_link_libraries(testBlockIntegration CoreLib UtilsLib)

add_executable(testBlockUnit Core/unit/test_block_unit.cpp)
target_link_libraries(testBlockUnit CoreLib UtilsLib)

add_executable(testBlockAccessorIntegration Core/integration/test_blockaccessor_integration.cpp)
target_link_libraries(testBlockAccessorIntegration CoreLib UtilsLib)

add_executable(testBlockAccessorUnit Core/unit/test_blockaccessor_unit.cpp)
target_link_libraries(testBlockAccessorUnit CoreLib UtilsLib)

add_executable(testCorrelationRepoIntegration Core/integration/test_correlationrepo_integration.cpp)
target_link_libraries(testCorrelationRepoIntegration CoreLib UtilsLib)

add_executable(testCorrelationRepoUnit Core/unit/test_correlationrepo_unit.cpp)
target_link_libraries(testCorrelationRepoUnit CoreLib UtilsLib)

add_executable(testParameterRouterIntegration Core/integration/test_parameterrouter_integration.cpp)
target_link_libraries(testParameterRouterIntegration CoreLib UtilsLib)

add_executable(testParameterRouterUnit Core/unit/test_parameterrouter_unit.cpp)
target_link_libraries(testParameterRouterUnit CoreLib UtilsLib)

add_executable(testDBMementoIntegration Core/integration/test_dbmemento_integration.cpp)
target_link_libraries(testDBMementoIntegration CoreLib UtilsLib)

add_executable(testDBMementoUnit Core/unit/test_dbmemento_unit.cpp)
target_link_libraries(testDBMementoUnit CoreLib UtilsLib)

add_executable(testMemoryManagerIntegration Core/integration/test_memory_manager_integration.cpp)
target_link_libraries(testMemoryManagerIntegration CoreLib UtilsLib)

add_executable(testMemoryManagerUnit Core/unit/test_memory_manager_unit.cpp)
target_link_libraries(testMemoryManagerUnit CoreLib UtilsLib)

target_include_directories(testMemoryManagerIntegration PRIVATE Core/support)
target_include_directories(testMemoryManagerUnit PRIVATE Core/support)

add_executable(testParametersIntegration Core/integration/test_parameters_integration.cpp)
target_link_libraries(testParametersIntegration CoreLib UtilsLib)

add_executable(testParametersUnit Core/unit/test_parameters_unit.cpp)
target_link_libraries(testParametersUnit CoreLib UtilsLib)

target_include_directories(testParametersIntegration PRIVATE Core/support)
target_include_directories(testParametersUnit PRIVATE Core/support)

add_executable(testSMparams Core/testSMParameters.cpp)
target_link_libraries(testSMparams CoreLib)

#MartyInterface

add_executable(testCompilerStrategyUnit ExternalIntegration/MartyInterface/unit/test_compiler_strategy_unit.cpp)
target_link_libraries(testCompilerStrategyUnit MartyInterfaceLib UtilsLib)

add_executable(testCompilerStrategyIntegration ExternalIntegration/MartyInterface/integration/test_compiler_strategy_integration.cpp)
target_link_libraries(testCompilerStrategyIntegration MartyInterfaceLib UtilsLib)

add_executable(testModelFileChecker ExternalIntegration/MartyInterface/unit/test_model_file_checker_unit.cpp)
target_link_libraries(testModelFileChecker MartyInterfaceLib UtilsLib)

add_executable(testGeneralModelModifierUnit ExternalIntegration/MartyInterface/unit/test_general_model_modifier_unit.cpp)
target_link_libraries(testGeneralModelModifierUnit MartyInterfaceLib UtilsLib)

add_executable(testGeneralModelModifierIntegration ExternalIntegration/MartyInterface/integration/test_general_model_modifier_integration.cpp)
target_link_libraries(testGeneralModelModifierIntegration MartyInterfaceLib UtilsLib)

add_executable(testExtractorUnit ExternalIntegration/MartyInterface/unit/test_extractor_unit.cpp)
target_link_libraries(testExtractorUnit MartyInterfaceLib UtilsLib)

add_executable(testFileWriterUnit ExternalIntegration/MartyInterface/unit/test_file_writer_unit.cpp)
target_link_libraries(testFileWriterUnit MartyInterfaceLib UtilsLib)

add_executable(testLineProcessorUnit ExternalIntegration/MartyInterface/unit/test_line_processor_unit.cpp)
target_link_libraries(testLineProcessorUnit MartyInterfaceLib UtilsLib)

add_executable(testParamWriterUnit ExternalIntegration/MartyInterface/unit/test_param_writer_unit.cpp)
target_link_libraries(testParamWriterUnit MartyInterfaceLib UtilsLib)

add_executable(testSMParamSetterUnit ExternalIntegration/MartyInterface/unit/test_param_setter_unit.cpp)
target_link_libraries(testSMParamSetterUnit MartyInterfaceLib UtilsLib)

add_executable(testWriterPipelineIntegration ExternalIntegration/MartyInterface/integration/test_writer_pipeline_integration.cpp)
target_link_libraries(testWriterPipelineIntegration MartyInterfaceLib UtilsLib)

#PhysicalModel

add_executable(testWilsonUnit PhysicalModel/unit/test_wilson_unit.cpp)
target_link_libraries(testWilsonUnit PhysicalModelLib UtilsLib)

add_executable(testWilsonIntegration PhysicalModel/integration/test_wilson_integration.cpp)
target_link_libraries(testWilsonIntegration PhysicalModelLib UtilsLib)

add_executable(testMartyWilsonUnit PhysicalModel/unit/test_marty_wilson_unit.cpp)
target_link_libraries(testMartyWilsonUnit PhysicalModelLib UtilsLib)

add_executable(testMartyWilsonIntegration PhysicalModel/integration/test_marty_wilson_integration.cpp)
target_link_libraries(testMartyWilsonIntegration PhysicalModelLib UtilsLib)

add_executable(testWilsonGroupUnit PhysicalModel/unit/test_wilsongroup_unit.cpp)
target_link_libraries(testWilsonGroupUnit PhysicalModelLib UtilsLib)

add_executable(testWilsonGroupIntegration PhysicalModel/integration/test_wilsongroup_integration.cpp)
target_link_libraries(testWilsonGroupIntegration PhysicalModelLib UtilsLib)

add_executable(testWilsonRegistryUnit PhysicalModel/unit/test_wilson_registry_unit.cpp)
target_link_libraries(testWilsonRegistryUnit PhysicalModelLib UtilsLib)

add_executable(testWilsonRegistryIntegration PhysicalModel/integration/test_wilson_registry_integration.cpp)
target_link_libraries(testWilsonRegistryIntegration PhysicalModelLib UtilsLib)

add_executable(testCoefficientGroupBuilderUnit PhysicalModel/unit/test_coefficient_group_builder_unit.cpp)
target_link_libraries(testCoefficientGroupBuilderUnit PhysicalModelLib UtilsLib)

add_executable(testCoefficientGroupBuilderIntegration PhysicalModel/integration/test_coefficient_group_builder_integration.cpp)
target_link_libraries(testCoefficientGroupBuilderIntegration PhysicalModelLib UtilsLib)

add_executable(testGroupDefinitionUnit PhysicalModel/unit/test_group_definition_unit.cpp)
target_link_libraries(testGroupDefinitionUnit PhysicalModelLib UtilsLib)

add_executable(testGroupDefinitionIntegration PhysicalModel/integration/test_group_definition_integration.cpp)
target_link_libraries(testGroupDefinitionIntegration PhysicalModelLib UtilsLib)

add_executable(testCustomWilsonUnit PhysicalModel/unit/test_custom_wilson_unit.cpp)
target_link_libraries(testCustomWilsonUnit PhysicalModelLib UtilsLib)

add_executable(testCustomWilsonIntegration PhysicalModel/integration/test_custom_wilson_integration.cpp)
target_link_libraries(testCustomWilsonIntegration PhysicalModelLib UtilsLib)

add_executable(testWilsonManagerUnit PhysicalModel/unit/test_wilson_manager_unit.cpp)
target_link_libraries(testWilsonManagerUnit PhysicalModelLib UtilsLib)

add_executable(testWilsonManagerIntegration PhysicalModel/integration/test_wilson_manager_integration.cpp)
target_link_libraries(testWilsonManagerIntegration PhysicalModelLib UtilsLib)

add_executable(testHasWilsonUnit PhysicalModel/unit/test_has_wilson_unit.cpp)
target_link_libraries(testHasWilsonUnit PhysicalModelLib UtilsLib)

add_executable(testHasWilsonIntegration PhysicalModel/integration/test_has_wilson_integration.cpp)
target_link_libraries(testHasWilsonIntegration PhysicalModelLib UtilsLib)

#BusinessLogic

add_executable(testObsWilsonHelperUnit BusinessLogic/unit/test_obswilsonhelper_unit.cpp)
target_link_libraries(testObsWilsonHelperUnit BusinessLogicLib UtilsLib)

add_executable(testObsWilsonHelperIntegration BusinessLogic/integration/test_obswilsonhelper_integration.cpp)
target_link_libraries(testObsWilsonHelperIntegration BusinessLogicLib UtilsLib)

add_executable(testDecayParentUnit BusinessLogic/unit/test_decayparent_unit.cpp)
target_link_libraries(testDecayParentUnit BusinessLogicLib UtilsLib)

add_executable(testDecayParentIntegration BusinessLogic/integration/test_decayparent_integration.cpp)
target_link_libraries(testDecayParentIntegration BusinessLogicLib UtilsLib)

add_executable(testObservableUnit BusinessLogic/unit/test_observable_unit.cpp)
target_link_libraries(testObservableUnit BusinessLogicLib UtilsLib)

add_executable(testObservableIntegration BusinessLogic/integration/test_observable_integration.cpp)
target_link_libraries(testObservableIntegration BusinessLogicLib UtilsLib)

add_executable(testObsManagerUnit BusinessLogic/unit/test_obsmanager_unit.cpp)
target_link_libraries(testObsManagerUnit BusinessLogicLib UtilsLib)

add_executable(testObsManagerIntegration BusinessLogic/integration/test_obsmanager_integration.cpp)
target_link_libraries(testObsManagerIntegration BusinessLogicLib UtilsLib)

#Statistic

add_executable(testMarginalsUnit Statistic/unit/test_marginals_unit.cpp)
target_link_libraries(testMarginalsUnit StatisticLib UtilsLib)

add_executable(testMarginalsIntegration Statistic/integration/test_marginals_integration.cpp)
target_link_libraries(testMarginalsIntegration StatisticLib UtilsLib)

add_executable(testCopulasUnit Statistic/unit/test_copulas_unit.cpp)
target_link_libraries(testCopulasUnit StatisticLib UtilsLib)

add_executable(testCopulasIntegration Statistic/integration/test_copulas_integration.cpp)
target_link_libraries(testCopulasIntegration StatisticLib UtilsLib)

add_executable(testJointDistributionUnit Statistic/unit/test_jointdistribution_unit.cpp)
target_link_libraries(testJointDistributionUnit StatisticLib UtilsLib)

add_executable(testJointDistributionIntegration Statistic/integration/test_jointdistribution_integration.cpp)
target_link_libraries(testJointDistributionIntegration StatisticLib UtilsLib)

add_executable(testLikelihoodBaseUnit Statistic/unit/test_likelihood_base_unit.cpp)
target_link_libraries(testLikelihoodBaseUnit StatisticLib UtilsLib)

add_executable(testLikelihoodBaseIntegration Statistic/integration/test_likelihood_base_integration.cpp)
target_link_libraries(testLikelihoodBaseIntegration StatisticLib UtilsLib)

add_executable(testGradientHelperUnit Statistic/unit/test_gradienthelper_unit.cpp)
target_link_libraries(testGradientHelperUnit StatisticLib UtilsLib)

add_executable(testGradientHelperIntegration Statistic/integration/test_gradienthelper_integration.cpp)
target_link_libraries(testGradientHelperIntegration StatisticLib UtilsLib)

add_executable(testProfilerUnit Statistic/unit/test_profiler_unit.cpp)
target_link_libraries(testProfilerUnit StatisticLib UtilsLib)

add_executable(testProfilerIntegration Statistic/integration/test_profiler_integration.cpp)
target_link_libraries(testProfilerIntegration StatisticLib UtilsLib)

add_executable(testIndexingUnit Statistic/unit/test_indexing_unit.cpp)
target_link_libraries(testIndexingUnit StatisticLib UtilsLib)

add_executable(testIndexingIntegration Statistic/integration/test_indexing_integration.cpp)
target_link_libraries(testIndexingIntegration StatisticLib UtilsLib)

add_executable(testChiSquaredLikelihoodUnit Statistic/unit/test_chisquaredlikelihood_unit.cpp)
target_link_libraries(testChiSquaredLikelihoodUnit StatisticLib UtilsLib)

add_executable(testChiSquaredLikelihoodIntegration Statistic/integration/test_chisquaredlikelihood_integration.cpp)
target_link_libraries(testChiSquaredLikelihoodIntegration StatisticLib UtilsLib)

add_executable(testProfilingStrategyUnit Statistic/unit/test_profilingstrategy_unit.cpp)
target_link_libraries(testProfilingStrategyUnit StatisticLib UtilsLib)

add_executable(testProfilingStrategyIntegration Statistic/integration/test_profilingstrategy_integration.cpp)
target_link_libraries(testProfilingStrategyIntegration StatisticLib UtilsLib)

add_executable(testProfiledLikelihood2DUnit Statistic/unit/test_profiledlikelihood2d_unit.cpp)
target_link_libraries(testProfiledLikelihood2DUnit StatisticLib UtilsLib)

add_executable(testProfiledLikelihood2DIntegration Statistic/integration/test_profiledlikelihood2d_integration.cpp)
target_link_libraries(testProfiledLikelihood2DIntegration StatisticLib UtilsLib)

add_executable(testWithGaussianConstraintsUnit Statistic/unit/test_withgaussianconstraints_unit.cpp)
target_link_libraries(testWithGaussianConstraintsUnit StatisticLib UtilsLib)

add_executable(testWithGaussianConstraintsIntegration Statistic/integration/test_withgaussianconstraints_integration.cpp)
target_link_libraries(testWithGaussianConstraintsIntegration StatisticLib UtilsLib)
# --------- Wilson ---------


add_executable(testSMWilson Wilson/testSMWilson.cpp)
target_link_libraries(testSMWilson PhysicalModelLib UtilsLib)

add_executable(testSMCQPrimeWilson Wilson/testSMCQPrimeWilson.cpp)
target_link_libraries(testSMCQPrimeWilson PhysicalModelLib UtilsLib)

add_executable(testBase1Wilson Wilson/testBase1Wilson.cpp)
target_link_libraries(testBase1Wilson PhysicalModelLib UtilsLib)

add_executable(testBase2Wilson Wilson/testBase2Wilson.cpp)
target_link_libraries(testBase2Wilson PhysicalModelLib UtilsLib)

if(BUILD_WITH_2HDMC)
    add_executable(testTHDMWilson Wilson/testTHDMWilson.cpp)
    target_link_libraries(testTHDMWilson PhysicalModelLib UtilsLib)

    add_executable(testTHDMCQPrimeWilson Wilson/testTHDMCQPrimeWilson.cpp)
    target_link_libraries(testTHDMCQPrimeWilson PhysicalModelLib UtilsLib)
endif()

if(BUILD_WITH_SOFTSUSY)
    add_executable(testSUSYWilson Wilson/testSUSYWilson.cpp)
    target_link_libraries(testSUSYWilson PhysicalModelLib UtilsLib)

    add_executable(testSUSYCQPrimeWilson Wilson/testSUSYCQPrimeWilson.cpp)
    target_link_libraries(testSUSYCQPrimeWilson PhysicalModelLib UtilsLib)
endif()

#Observables

add_executable(testObservables Observables/test_observables.cpp)
target_link_libraries(testObservables BusinessLogicLib UtilsLib)

include(CTest)


add_test(NAME TestLogger COMMAND testCommon)
set_tests_properties(TestLogger PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

#Common

add_test(NAME testCommon COMMAND testCommon)
set_tests_properties(testCommon PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testlhaidunit COMMAND testlhaidunit)
set_tests_properties(testlhaidunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testlhaidintegration COMMAND testlhaidintegration)
set_tests_properties(testlhaidintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testblocknameunit COMMAND testblocknameunit)
set_tests_properties(testblocknameunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testblocknameintegration COMMAND testblocknameintegration)
set_tests_properties(testblocknameintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testparamidunit COMMAND testparamidunit)
set_tests_properties(testparamidunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testparamidintegration COMMAND testparamidintegration)
set_tests_properties(testparamidintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testenummapperunit COMMAND testenummapperunit)
set_tests_properties(testenummapperunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

add_test(NAME testenummapperintegration COMMAND testenummapperintegration)
set_tests_properties(testenummapperintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")

#DataBase

add_test(NAME testReader COMMAND testReader)
set_tests_properties(testReader PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testjsonreader COMMAND testjsonreader)
set_tests_properties(testjsonreader PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testmappingjson COMMAND testmappingjson)
set_tests_properties(testmappingjson PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME teststatjson COMMAND teststatjson)
set_tests_properties(teststatjson PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testlhaelementunit COMMAND testlhaelementunit)
set_tests_properties(testlhaelementunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testlhaelementintegration COMMAND testlhaelementintegration)
set_tests_properties(testlhaelementintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testdbnodeunit COMMAND testdbnodeunit)
set_tests_properties(testdbnodeunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testdbnodeintegration COMMAND testdbnodeintegration)
set_tests_properties(testdbnodeintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testdbmanagerunit COMMAND testdbmanagerunit)
set_tests_properties(testdbmanagerunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testdbmanagerintegration COMMAND testdbmanagerintegration)
set_tests_properties(testdbmanagerintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testlhablocksunit COMMAND testlhablocksunit)
set_tests_properties(testlhablocksunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testlhablocksintegration COMMAND testlhablocksintegration)
set_tests_properties(testlhablocksintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testjsonparserunit COMMAND testjsonparserunit)
set_tests_properties(testjsonparserunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testjsonparserintegration COMMAND testjsonparserintegration)
set_tests_properties(testjsonparserintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testlhaparserunit COMMAND testlhaparserunit)
set_tests_properties(testlhaparserunit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME testlhaparserintegration COMMAND testlhaparserintegration)
set_tests_properties(testlhaparserintegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

#Pandas

add_test(NAME testSerieUnit COMMAND testSerieUnit)
set_tests_properties(testSerieUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "pandas")

add_test(NAME testSerieIntegration COMMAND testSerieIntegration)
set_tests_properties(testSerieIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "pandas")

add_test(NAME testDataFrameUnit COMMAND testDataFrameUnit)
set_tests_properties(testDataFrameUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "pandas")

add_test(NAME testDataFrameIntegration COMMAND testDataFrameIntegration)
set_tests_properties(testDataFrameIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "pandas")

#Core

add_test(NAME testParameterIntegration COMMAND testParameterIntegration)
set_tests_properties(testParameterIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testParameterUnit COMMAND testParameterUnit)
set_tests_properties(testParameterUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testBlockIntegration COMMAND testBlockIntegration)
set_tests_properties(testBlockIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testBlockUnit COMMAND testBlockUnit)
set_tests_properties(testBlockUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testBlockAccessorIntegration COMMAND testBlockAccessorIntegration)
set_tests_properties(testBlockAccessorIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testBlockAccessorUnit COMMAND testBlockAccessorUnit)
set_tests_properties(testBlockAccessorUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testCorrelationRepoIntegration COMMAND testCorrelationRepoIntegration)
set_tests_properties(testCorrelationRepoIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testCorrelationRepoUnit COMMAND testCorrelationRepoUnit)
set_tests_properties(testCorrelationRepoUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testParameterRouterIntegration COMMAND testParameterRouterIntegration)
set_tests_properties(testParameterRouterIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testParameterRouterUnit COMMAND testParameterRouterUnit)
set_tests_properties(testParameterRouterUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testDBMementoIntegration COMMAND testDBMementoIntegration)
set_tests_properties(testDBMementoIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testDBMementoUnit COMMAND testDBMementoUnit)
set_tests_properties(testDBMementoUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testMemoryManagerIntegration COMMAND testMemoryManagerIntegration)
set_tests_properties(testMemoryManagerIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testMemoryManagerUnit COMMAND testMemoryManagerUnit)
set_tests_properties(testMemoryManagerUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testParametersIntegration COMMAND testParametersIntegration)
set_tests_properties(testParametersIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testParametersUnit COMMAND testParametersUnit)
set_tests_properties(testParametersUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

add_test(NAME testSMparams COMMAND testSMparams)
set_tests_properties(testSMparams PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")

#MartyInterface

add_test(NAME testCompilerStrategyUnit COMMAND testCompilerStrategyUnit)
set_tests_properties(testCompilerStrategyUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testCompilerStrategyIntegration COMMAND testCompilerStrategyIntegration)
set_tests_properties(testCompilerStrategyIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testModelFileChecker COMMAND testModelFileChecker)
set_tests_properties(testModelFileChecker PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testGeneralModelModifierUnit COMMAND testGeneralModelModifierUnit)
set_tests_properties(testGeneralModelModifierUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testGeneralModelModifierIntegration COMMAND testGeneralModelModifierIntegration)
set_tests_properties(testGeneralModelModifierIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testExtractorUnit COMMAND testExtractorUnit)
set_tests_properties(testExtractorUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testFileWriterUnit COMMAND testFileWriterUnit)
set_tests_properties(testFileWriterUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testLineProcessorUnit COMMAND testLineProcessorUnit)
set_tests_properties(testLineProcessorUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testParamWriterUnit COMMAND testParamWriterUnit)
set_tests_properties(testParamWriterUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testSMParamSetterUnit COMMAND testSMParamSetterUnit)
set_tests_properties(testSMParamSetterUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

add_test(NAME testWriterPipelineIntegration COMMAND testWriterPipelineIntegration)
set_tests_properties(testWriterPipelineIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "martyinterface")

#PhysicalModel

add_test(NAME testWilsonUnit COMMAND testWilsonUnit)
set_tests_properties(testWilsonUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testWilsonIntegration COMMAND testWilsonIntegration)
set_tests_properties(testWilsonIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testWilsonGroupUnit COMMAND testWilsonGroupUnit)
set_tests_properties(testWilsonGroupUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testWilsonGroupIntegration COMMAND testWilsonGroupIntegration)
set_tests_properties(testWilsonGroupIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")


add_test(NAME testWilsonRegistryUnit COMMAND testWilsonRegistryUnit)
set_tests_properties(testWilsonRegistryUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testWilsonRegistryIntegration COMMAND testWilsonRegistryIntegration)
set_tests_properties(testWilsonRegistryIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testCoefficientGroupBuilderUnit COMMAND testCoefficientGroupBuilderUnit)
set_tests_properties(testCoefficientGroupBuilderUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testCoefficientGroupBuilderIntegration COMMAND testCoefficientGroupBuilderIntegration)
set_tests_properties(testCoefficientGroupBuilderIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testGroupDefinitionUnit COMMAND testGroupDefinitionUnit)
set_tests_properties(testGroupDefinitionUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testGroupDefinitionIntegration COMMAND testGroupDefinitionIntegration)
set_tests_properties(testGroupDefinitionIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testCustomWilsonUnit COMMAND testCustomWilsonUnit)
set_tests_properties(testCustomWilsonUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testCustomWilsonIntegration COMMAND testCustomWilsonIntegration)
set_tests_properties(testCustomWilsonIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testWilsonManagerUnit COMMAND testWilsonManagerUnit)
set_tests_properties(testWilsonManagerUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testWilsonManagerIntegration COMMAND testWilsonManagerIntegration)
set_tests_properties(testWilsonManagerIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testHasWilsonUnit COMMAND testHasWilsonUnit)
set_tests_properties(testHasWilsonUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

add_test(NAME testHasWilsonIntegration COMMAND testHasWilsonIntegration)
set_tests_properties(testHasWilsonIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "physicalmodel")

#BusinessLogic

add_test(NAME testObsWilsonHelperUnit COMMAND testObsWilsonHelperUnit)
set_tests_properties(testObsWilsonHelperUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testObsWilsonHelperIntegration COMMAND testObsWilsonHelperIntegration)
set_tests_properties(testObsWilsonHelperIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testDecayParentUnit COMMAND testDecayParentUnit)
set_tests_properties(testDecayParentUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testDecayParentIntegration COMMAND testDecayParentIntegration)
set_tests_properties(testDecayParentIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testObservableUnit COMMAND testObservableUnit)
set_tests_properties(testObservableUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testObservableIntegration COMMAND testObservableIntegration)
set_tests_properties(testObservableIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testObsManagerUnit COMMAND testObsManagerUnit)
set_tests_properties(testObsManagerUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

add_test(NAME testObsManagerIntegration COMMAND testObsManagerIntegration)
set_tests_properties(testObsManagerIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "businesslogic")

#Statistic

add_test(NAME testMarginalsUnit COMMAND testMarginalsUnit)
set_tests_properties(testMarginalsUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testMarginalsIntegration COMMAND testMarginalsIntegration)
set_tests_properties(testMarginalsIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testCopulasUnit COMMAND testCopulasUnit)
set_tests_properties(testCopulasUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testCopulasIntegration COMMAND testCopulasIntegration)
set_tests_properties(testCopulasIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testJointDistributionUnit COMMAND testJointDistributionUnit)
set_tests_properties(testJointDistributionUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testJointDistributionIntegration COMMAND testJointDistributionIntegration)
set_tests_properties(testJointDistributionIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testLikelihoodBaseUnit COMMAND testLikelihoodBaseUnit)
set_tests_properties(testLikelihoodBaseUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testLikelihoodBaseIntegration COMMAND testLikelihoodBaseIntegration)
set_tests_properties(testLikelihoodBaseIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testGradientHelperUnit COMMAND testGradientHelperUnit)
set_tests_properties(testGradientHelperUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testGradientHelperIntegration COMMAND testGradientHelperIntegration)
set_tests_properties(testGradientHelperIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testProfilerUnit COMMAND testProfilerUnit)
set_tests_properties(testProfilerUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testProfilerIntegration COMMAND testProfilerIntegration)
set_tests_properties(testProfilerIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testIndexingUnit COMMAND testIndexingUnit)
set_tests_properties(testIndexingUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testIndexingIntegration COMMAND testIndexingIntegration)
set_tests_properties(testIndexingIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testChiSquaredLikelihoodUnit COMMAND testChiSquaredLikelihoodUnit)
set_tests_properties(testChiSquaredLikelihoodUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testChiSquaredLikelihoodIntegration COMMAND testChiSquaredLikelihoodIntegration)
set_tests_properties(testChiSquaredLikelihoodIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testProfilingStrategyUnit COMMAND testProfilingStrategyUnit)
set_tests_properties(testProfilingStrategyUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testProfilingStrategyIntegration COMMAND testProfilingStrategyIntegration)
set_tests_properties(testProfilingStrategyIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testProfiledLikelihood2DUnit COMMAND testProfiledLikelihood2DUnit)
set_tests_properties(testProfiledLikelihood2DUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testProfiledLikelihood2DIntegration COMMAND testProfiledLikelihood2DIntegration)
set_tests_properties(testProfiledLikelihood2DIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testWithGaussianConstraintsUnit COMMAND testWithGaussianConstraintsUnit)
set_tests_properties(testWithGaussianConstraintsUnit PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

add_test(NAME testWithGaussianConstraintsIntegration COMMAND testWithGaussianConstraintsIntegration)
set_tests_properties(testWithGaussianConstraintsIntegration PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "statistic")

#other testings
add_test(NAME TestLHAReader COMMAND testReader)
set_tests_properties(TestLHAReader PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME TestJsonReader COMMAND testjsonreader)
set_tests_properties(TestJsonReader PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME TestMappingJson COMMAND testmappingjson)
set_tests_properties(TestMappingJson PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")

add_test(NAME TestStatJson COMMAND teststatjson)
set_tests_properties(TestStatJson PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "database")


add_test(NAME TestSMParameters COMMAND testSMparams)
set_tests_properties(TestSMParameters PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "core")


add_test(NAME WCTestSM COMMAND testSMWilson)
set_tests_properties(WCTestSM PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")

add_test(NAME WCTestSMCQPrime COMMAND testSMCQPrimeWilson)
set_tests_properties(WCTestSMCQPrime PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")

add_test(NAME WCTestBase1 COMMAND testBase1Wilson)
set_tests_properties(WCTestBase1 PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")

add_test(NAME WCTestBase2 COMMAND testBase2Wilson)
set_tests_properties(WCTestBase2 PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")

if(BUILD_WITH_2HDMC)
    add_test(NAME WCTestTHDM COMMAND testTHDMWilson)
    set_tests_properties(WCTestTHDM PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")

    add_test(NAME WCTestTHDMCQPrime COMMAND testTHDMCQPrimeWilson)
    set_tests_properties(WCTestTHDMCQPrime PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")
else()
    message("Couldn't do the THDM test, 2HDMC is not installed or not configured, use -DBUILD_WITH_2HDMC=ON")
endif()

if(BUILD_WITH_SOFTSUSY)
    add_test(NAME WCTestSUSY COMMAND testSUSYWilson)
    set_tests_properties(WCTestSUSY PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")

    add_test(NAME WCTestSUSYCQPrime COMMAND testSUSYCQPrimeWilson)
    set_tests_properties(WCTestSUSYCQPrime PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "wilson")
else()
    message("Couldn't do the SUSY test, SoftSusy is not installed or not configured, use -DBUILD_WITH_SOFTSUSY=ON")
endif()
