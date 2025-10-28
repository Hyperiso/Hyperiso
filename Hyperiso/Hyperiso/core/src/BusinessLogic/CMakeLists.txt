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
    domain/Calculators/BPFFCalculator.cpp
    domain/Calculators/BVFFCalculator.cpp
    domain/Calculators/BaseQCDfCalculator.cpp
    domain/Calculators/BPQCDfCalculator.cpp
    domain/Calculators/BVQCDfCalculator.cpp
    domain/Decays/BKsllDecay.cpp
    domain/Decays/BKstarGammaDecay.cpp
    domain/Decays/BllDecay.cpp
    domain/Decays/BXsDecay.cpp
    domain/Decays/BlnuDecay.cpp
    domain/Decays/BDlnuDecay.cpp
    domain/Decays/BDstarlnuDecay.cpp
    domain/Decays/BXsllDecay.cpp
    domain/Decays/M0_Mixing.cpp
    domain/Decays/BsPhiDecay.cpp
    domain/Decays/BKllDecay.cpp
    domain/Decays/KllDecay.cpp
    # domain/Wilson/WilsonObs.cpp
    domain/Compound.cpp
    domain/Observable.cpp
    domain/DecayParent.cpp
    domain/Node.cpp
    domain/ModelEvaluator.cpp
    domain/ObsManager.cpp
    domain/ObsWilsonHelper.cpp
    adapters/CorrelationProxy.cpp
    adapters/ObsParameterProxy.cpp
    adapters/ObsParameterMutator.cpp
    adapters/ObsQCDProxy.cpp
    adapters/ObsWilsonBuilder.cpp
    adapters/ObsWilsonProxy.cpp
    adapters/WilsonFreezer.cpp
    adapters/ObservableInterface.cpp
)
set_target_properties(BusinessLogicLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# set_target_properties(BusinessLogicLib PROPERTIES PUBLIC_HEADER "Decays.h;Compound.h;DecayParent.h;Node.h;Observable.h;ModelEvaluator.h;ObsManager.h;ObservableInterface.h;")

# target_include_directories(BusinessLogicLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/domain)
# target_include_directories(BusinessLogicLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ports)
# target_include_directories(BusinessLogicLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
# target_include_directories(BusinessLogicLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/domain/Decays)

target_include_directories(BusinessLogicLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Calculators>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Decays>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic/domain/Calculators>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic/domain/Decays>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic/ports>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic/adapters>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(BusinessLogicLib PUBLIC PhysicalModelLib CoreLib CommonLib ConfiguredModule MathLib)

add_executable(TEST app/main.cpp)
target_link_libraries(TEST BusinessLogicLib)

add_executable(CAMILIA app/main_camilia.cpp)
target_link_libraries(CAMILIA BusinessLogicLib)

add_executable(BSLL_FF app/test_bsll_ff.cpp)
target_link_libraries(BSLL_FF BusinessLogicLib)

add_executable(TEST_CUSTOM_OBS app/main_custom_obs.cpp)
target_link_libraries(TEST_CUSTOM_OBS BusinessLogicLib)

include(GNUInstallDirs)

install(TARGETS BusinessLogicLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/BusinessLogic
  FILES_MATCHING
    PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)