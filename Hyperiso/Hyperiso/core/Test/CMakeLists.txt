cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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

add_executable(testWilsonGroupUnit PhysicalModel/unit/test_wilsongroup_unit.cpp)
target_link_libraries(testWilsonGroupUnit PhysicalModelLib UtilsLib)

add_executable(testWilsonGroupIntegration PhysicalModel/integration/test_wilsongroup_integration.cpp)
target_link_libraries(testWilsonGroupIntegration PhysicalModelLib UtilsLib)

#Wilson


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



