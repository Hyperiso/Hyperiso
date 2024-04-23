cmake_minimum_required(VERSION 3.20)
project(TestHyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CONFIG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../config.json")
get_filename_component(PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../" ABSOLUTE)
set(CONFIG_CONTENT "{\"project_root\": \"${PROJECT_ROOT}\"}")
file(WRITE "${CONFIG_FILE}" "${CONFIG_CONTENT}")
message(STATUS "Configuration file created : ${CONFIG_FILE}")

# Option for building with SOFTSUSY
option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)



# Define common source files
file(GLOB COMMON_SOURCES
     "../Core/*.cpp"
     "../DataBase/*.cpp"
     # "../logique_metier/*.cpp"
     "../Math/*.cpp"
     "../Physical_Model/*.cpp"
     "../External_Integration/SoftSusy.cpp"
     )
file(GLOB SUSY "../Physical_Model/BSM/*.cpp" "../Math/Math_BSM/*.cpp")
set(TEST_READER_SOURCE "testReader.cpp")
set(TEST_WILSON_SOURCE "testWilson.cpp")
set(TEST_SUSYWILSON_SOURCE "testSusyWilson.cpp" )
set(TEST_PARAMS_SOURCE "testParams.cpp")
set(TEST_QCDPARAMETERS_SOURCE "testQCDParameters.cpp") # Ajout de la nouvelle source

add_executable(testReader ${COMMON_SOURCES} ${TEST_READER_SOURCE})
add_executable(testWilson ${COMMON_SOURCES} ${TEST_WILSON_SOURCE})
add_executable(testSusyWilson ${COMMON_SOURCES} ${SUSY} ${TEST_SUSYWILSON_SOURCE})
add_executable(testParams ${COMMON_SOURCES} ${TEST_PARAMS_SOURCE})
add_executable(testQCDParameters ${COMMON_SOURCES} ${TEST_QCDPARAMETERS_SOURCE}) # Ajout de la nouvelle cible

# Add include directories for testReader
target_include_directories(testReader 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    # "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
    "${CMAKE_SOURCE_DIR}/../External_Integration"
)

# Add include directories for testWilson
target_include_directories(testWilson 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    # "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
    "${CMAKE_SOURCE_DIR}/../External_Integration"
)

# Add include directories for testWilson
target_include_directories(testParams
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
    "${CMAKE_SOURCE_DIR}/../External_Integration"
)

# Add include directories for testSusyWilson
target_include_directories(testSusyWilson
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Math/Math_BSM"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
    "${CMAKE_SOURCE_DIR}/../Physical_Model/BSM"
    "${CMAKE_SOURCE_DIR}/../External_Integration"
)

# Add include directories for testQCDParameters
target_include_directories(testQCDParameters 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/../Core"
    "${CMAKE_SOURCE_DIR}/../DataBase"
    # "${CMAKE_SOURCE_DIR}/../logique_metier"
    "${CMAKE_SOURCE_DIR}/../Math"
    "${CMAKE_SOURCE_DIR}/../Physical_Model"
    "${CMAKE_SOURCE_DIR}/../External_Integration"
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

if(BUILD_WITH_SOFTSUSY)
  # Include the ExternalProject module
  include(ExternalProject)
  
  # Add SOFTSUSY as an external project
  ExternalProject_Add(SOFTSUSY
    GIT_REPOSITORY "https://github.com/BAllanach/softsusy.git" # Specify the git repository of SOFTSUSY
    GIT_TAG "v4.1.19" # Specify a commit, tag, or branch
    PREFIX "${CMAKE_BINARY_DIR}/../../External_Integration/SOFTSUSY"
    # CONFIGURE_COMMAND "autoreconf -fiv && ./configure"
    CONFIGURE_COMMAND ""
    COMMAND autoreconf -fiv
    COMMAND ./configure
    BUILD_COMMAND "make"
    BUILD_IN_SOURCE TRUE
    INSTALL_COMMAND ""
  )
  
  # Specify SOFTSUSY's build directory as an include directory
  # This assumes SOFTSUSY generates header files or you need to include something from its source.
  # Adjust the path according to where the needed headers are located.
#   externalproject_get_property(SOFTSUSY SOURCE_DIR)
#   target_include_directories(test PUBLIC "${SOURCE_DIR}")
  ExternalProject_Get_Property(SOFTSUSY SOURCE_DIR)
  target_include_directories(testParams PUBLIC "${SOURCE_DIR}")
  add_dependencies(testParams SOFTSUSY)
endif()