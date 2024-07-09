cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT TARGET UtilsLib)
    add_subdirectory(Utils UtilsLib)
endif()

include(Catch)
set(CTEST_ENVIRONMENT_VARIABLE "CTEST_RUN")
# add_executable(testSMWilson testSMWilson.cpp)
# target_link_libraries(testSMWilson PhysicalModelLib)

# add_executable(testSusyWilson testSusyWilson.cpp)
# target_link_libraries(testSusyWilson PhysicalModelLib)

# add_executable(testTHDMWilson testTHDMWilson.cpp)
# target_link_libraries(testTHDMWilson PhysicalModelLib)

add_executable(testWilson testWilson.cpp)
target_link_libraries(testWilson PhysicalModelLib)

add_executable(testObservable testObservable.cpp)
target_link_libraries(testObservable BusinessLogicLib)

add_executable(testParams testParams.cpp)
target_link_libraries(testParams CoreLib)

add_executable(testPrimeCQWilson testPrimeCQWilson.cpp)
target_link_libraries(testPrimeCQWilson PhysicalModelLib)

add_executable(testPrimeCQTHDMWilson testPrimeCQTHDMWilson.cpp)
target_link_libraries(testPrimeCQTHDMWilson PhysicalModelLib)

add_executable(testPrimeCQSUSYWilson testPrimeCQSUSYWilson.cpp)
target_link_libraries(testPrimeCQSUSYWilson PhysicalModelLib)
# add_executable(test2HDMC test2HDMC.cpp)
# target_link_libraries(test2HDMC ExternalIntegrationLib)

# add_executable(testQCDParameters testQCDParameters.cpp)
# target_link_libraries(testQCDParameters CoreLib)

add_executable(testReader testReader.cpp)
target_link_libraries(testReader DataBaseLib)

add_executable(testSMWilson Wilson/testSMWilson.cpp)
target_link_libraries(testSMWilson PhysicalModelLib UtilsLib)

add_executable(testSMCQPrimeWilson Wilson/testSMCQPrimeWilson.cpp)
target_link_libraries(testSMCQPrimeWilson PhysicalModelLib UtilsLib)

add_executable(testTHDMWilson Wilson/testTHDMWilson.cpp)
target_link_libraries(testTHDMWilson PhysicalModelLib UtilsLib)

add_executable(testTHDMCQPrimeWilson Wilson/testTHDMCQPrimeWilson.cpp)
target_link_libraries(testTHDMCQPrimeWilson PhysicalModelLib UtilsLib)

add_executable(testSUSYWilson Wilson/testSUSYWilson.cpp)
target_link_libraries(testSUSYWilson PhysicalModelLib UtilsLib)

add_executable(testSUSYCQPrimeWilson Wilson/testSUSYCQPrimeWilson.cpp)
target_link_libraries(testSUSYCQPrimeWilson PhysicalModelLib UtilsLib)


include(CTest)


add_test(NAME WCTestSM COMMAND testSMWilson)
set_tests_properties(WCTestSM PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1")

add_test(NAME WCTestSMCQPrime COMMAND testSMCQPrimeWilson)
set_tests_properties(WCTestSMCQPrime PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1")

add_test(NAME WCTestTHDM COMMAND testTHDMWilson)
set_tests_properties(WCTestTHDM PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1")

add_test(NAME WCTestTHDMCQPrime COMMAND testTHDMCQPrimeWilson)
set_tests_properties(WCTestTHDMCQPrime PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1")


add_test(NAME WCTestSUSY COMMAND testSUSYWilson)
set_tests_properties(WCTestSUSY PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1")

add_test(NAME WCTestSUSYCQPrime COMMAND testSUSYCQPrimeWilson)
set_tests_properties(WCTestSUSYCQPrime PROPERTIES ENVIRONMENT "DISABLE_LOGGER=1")
