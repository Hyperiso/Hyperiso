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

include(Catch)
set(CTEST_ENVIRONMENT_VARIABLE "CTEST_RUN")
# add_executable(testSMWilson testSMWilson.cpp)
# target_link_libraries(testSMWilson PhysicalModelLib)

# add_executable(testSusyWilson testSusyWilson.cpp)
# target_link_libraries(testSusyWilson PhysicalModelLib)

# add_executable(testTHDMWilson testTHDMWilson.cpp)
# target_link_libraries(testTHDMWilson PhysicalModelLib)

# add_executable(testWilson testWilson.cpp)
# target_link_libraries(testWilson PhysicalModelLib)

add_executable(testObservable testObservable.cpp)
target_link_libraries(testObservable BusinessLogicLib)

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

add_executable(testCommon Common/testLogger.cpp)
target_link_libraries(testCommon CommonLib)


add_executable(testReader DataBase/testReader.cpp)
target_link_libraries(testReader DataBaseLib ConfiguredModule)

add_executable(testjsonreader DataBase/testJsonReader.cpp)
target_link_libraries(testjsonreader DataBaseLib ConfiguredModule)


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



