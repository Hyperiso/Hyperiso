cmake_minimum_required(VERSION 3.20)
project(Hyperiso)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Option for building with SOFTSUSY
option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)

if(BUILD_WITH_SOFTSUSY)
  # Include the ExternalProject module
  include(ExternalProject)
  
  # Add SOFTSUSY as an external project
  ExternalProject_Add(SOFTSUSY
    GIT_REPOSITORY "https://github.com/BAllanach/softsusy.git" # Specify the git repository of SOFTSUSY
    GIT_TAG "v4.1.19" # Specify a commit, tag, or branch
    PREFIX "${CMAKE_BINARY_DIR}/External_Integration/SOFTSUSY"
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
endif()


file(GLOB SOURCES "./Core/*.cpp" "./DataBase/*.cpp" "./logique_metier/*.cpp" "./Math/*.cpp" "./Physical_Model/*.cpp" "./External_Integration/main_softsusy.cpp" "./External_Integration/SoftSusy.cpp")

add_executable(test ${SOURCES})

# Add include directories
target_include_directories(test 
    PUBLIC 
    "${CMAKE_SOURCE_DIR}/Core"
    "${CMAKE_SOURCE_DIR}/DataBase"
    "${CMAKE_SOURCE_DIR}/logique_metier"
    "${CMAKE_SOURCE_DIR}/Math"
    "${CMAKE_SOURCE_DIR}/Physical_Model"
    "${CMAKE_SOURCE_DIR}/External_Integration"
)

# Link any necessary libraries to the test executable
# For example:
# target_link_libraries(tests PRIVATE gtest)

# Optionally, add a custom target for running tests with CTest
# enable_testing()
# add_test(NAME AllTests COMMAND tests)

