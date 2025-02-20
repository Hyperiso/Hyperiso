cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(BusinessLogic)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

if (NOT TARGET PhysicalModelLib)
    add_subdirectory(../PhysicalModel PhysicalModelLib)
endif()
if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()
if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured ConfiguredModule)
endif()
if (NOT TARGET MathLib)
    add_subdirectory(../Math MathLib)
endif()
if (NOT TARGET CoreLib)
    add_subdirectory(../Core CoreLib)
endif()

add_library(BusinessLogicLib STATIC 
    Decays/BKstarDecay.cpp
    Decays/BllDecay.cpp
    Decays/BXsDecay.cpp
    Decays/BlnuDecay.cpp
    Decays/BDlnuDecay.cpp
    Decays/BDstarlnuDecay.cpp
    Compound.cpp
    DecayParent.cpp
    Node.cpp
    ModelEvaluator.cpp
    ObsManager.cpp
)
set_target_properties(BusinessLogicLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_target_properties(BusinessLogicLib PROPERTIES PUBLIC_HEADER "Decays.h;Compound.h;DecayParent.h;Node.h;Observable.h;ModelEvaluator.h;ObsManager.h;ObservableInterface.h;")

target_include_directories(BusinessLogicLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(BusinessLogicLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Decays)

target_link_libraries(BusinessLogicLib PUBLIC PhysicalModelLib CoreLib CommonLib ConfiguredModule MathLib)

add_executable(TEST main.cpp)
target_link_libraries(TEST BusinessLogicLib)