cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define common source files
file(GLOB COMMON_SOURCES
     "../Core/*.cpp"
     "../DataBase/*.cpp"
     # "../logique_metier/*.cpp"
     "../Math/*.cpp"
     "../Physical_Model/*.cpp")

set(TEST_READER_SOURCE "testReader.cpp")
set(TEST_WILSON_SOURCE "testWilson.cpp")
set(TEST_QCDPARAMETERS_SOURCE "testQCDParameters.cpp") # Ajout de la nouvelle source

add_executable(testReader ${COMMON_SOURCES} ${TEST_READER_SOURCE})
add_executable(testWilson ${COMMON_SOURCES} ${TEST_WILSON_SOURCE})
add_executable(testQCDParameters ${COMMON_SOURCES} ${TEST_QCDPARAMETERS_SOURCE}) # Ajout de la nouvelle cible

# Add include directories for testReader
target_include_directories(testReader 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    # "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
)

# Add include directories for testWilson
target_include_directories(testWilson 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    # "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
)

# Add include directories for testQCDParameters
target_include_directories(testQCDParameters 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    # "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
)

# Link any necessary libraries to the test executable
# For example:
# target_link_libraries(testReader PRIVATE gtest)
# target_link_libraries(testWilson PRIVATE gtest)
# target_link_libraries(testQCDParameters PRIVATE gtest) # Si nécessaire

# Optionally, add a custom target for running tests with CTest
# enable_testing()
# add_test(NAME ReaderTest COMMAND testReader)
# add_test(NAME WilsonTest COMMAND testWilson)
# add_test(NAME QCDParametersTest COMMAND testQCDParameters) # Ajout d'un test CTest

