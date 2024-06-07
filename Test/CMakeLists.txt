cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(testSMWilson testSMWilson.cpp)
target_link_libraries(testSMWilson PhysicalModelLib)

add_executable(testSusyWilson testSusyWilson.cpp)
target_link_libraries(testSusyWilson PhysicalModelLib)

add_executable(testTHDMWilson testTHDMWilson.cpp)
target_link_libraries(testTHDMWilson PhysicalModelLib)

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