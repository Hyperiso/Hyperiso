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

#Core


add_executable(testSMparams Core/testSMParameters.cpp)
target_link_libraries(testSMparams CoreLib)

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

include(CTest)


add_test(NAME TestLogger COMMAND testCommon)
set_tests_properties(TestLogger PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1" LABELS "common")


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



