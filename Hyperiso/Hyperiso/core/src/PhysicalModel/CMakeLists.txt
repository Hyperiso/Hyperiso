cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PhysicalModel)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

if (BUILD_WITH_2HDMC)
    message(STATUS "wow")
    add_compile_definitions(BUILD_WITH_2HDMC)
endif()

if (BUILD_WITH_SOFTSUSY)
    add_definitions(-BUILD_WITH_SOFTSUSY)
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

if (NOT TARGET MartyInterfaceLib)
    add_subdirectory(../ExternalIntegration/MartyInterface MartyInterfaceLib)
endif()


add_library(PhysicalModelLib STATIC 
    domain/MesonMixingRunningParameters.cpp
    domain/BWilsonRunningParameters.cpp
    domain/Wilson_parameters.cpp
    domain/Wilson.cpp
    domain/MartyWilson.cpp
    domain/BWilson.cpp
    domain/BWilsonGroup.cpp
    domain/ChargedCurrentWilson.cpp
    domain/ChargedCurrentsWilsonGroup.cpp
    domain/MesonMixingWilson.cpp
    domain/MesonMixingWilsonGroup.cpp
    domain/WilsonGroup.cpp
    domain/WilsonManager.cpp
    domain/MartyWilson.cpp
    adapters/ParameterProxy.cpp
    adapters/WilsonParamComposer.cpp
    adapters/ScaleSetter.cpp
    adapters/WilsonBuilder.cpp
    adapters/WilsonProvider.cpp
    adapters/BlockProxy.cpp
    BSM/epsilon_calculator.cpp
    BSM/susy_parameters.cpp
    BSM/thdm_parameters.cpp
    BSM/Wilson_SUSY_super.cpp
    BSM/Wilson_THDM_super.cpp
    )

set_target_properties(PhysicalModelLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(PhysicalModelLib PROPERTIES PUBLIC_HEADER "Interpolator.h;WilsonManager.h;Wilson.h;WilsonGroup.h;MartyWilson.h;Wilson_parameters.h;epsilon_calculator.h;susy_parameters.h;thdm_parameters;Wilson_susy.h;Wilson_THDM.h;WilsonManager.h")

# target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
# target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/ports)
# target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/domain)
# target_include_directories(PhysicalModelLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/BSM)


target_include_directories(PhysicalModelLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ports>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BSM>

    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/adapters>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/ports>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/BSM>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(PhysicalModelLib PUBLIC MartyInterfaceLib CoreLib CommonLib ConfiguredModule MathLib)

add_executable(test_wilson app/main_PM.cpp)
target_link_libraries(test_wilson PhysicalModelLib)

add_executable(test_wilson_marty app/main_marty.cpp)
target_link_libraries(test_wilson_marty PhysicalModelLib)

add_executable(test_contributions app/test_contributions.cpp)
target_link_libraries(test_contributions PhysicalModelLib)

add_executable(test_add_wilson app/main_phymod.cpp)
target_link_libraries(test_add_wilson PhysicalModelLib)

add_executable(camilia_test app/main_camilia.cpp)
target_link_libraries(camilia_test PhysicalModelLib)

add_executable(susy_test app/main_susy.cpp)
target_link_libraries(susy_test PhysicalModelLib)

add_executable(thdm_test app/main_thdm.cpp)
target_link_libraries(thdm_test PhysicalModelLib)

add_executable(TEST_SWITCH_LHA app/main_switch_lha.cpp)
target_link_libraries(TEST_SWITCH_LHA PhysicalModelLib)

add_executable(custom_test app/main_custom_wilson.cpp)
target_link_libraries(custom_test PhysicalModelLib)
# add_executable(test_martywilson main.cpp)
# target_link_libraries(test_martywilson PhysicalModelLib MartyInterfaceLib)

include(GNUInstallDirs)

install(TARGETS PhysicalModelLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel
        FILES_MATCHING
          PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
          PATTERN "*.cpp" EXCLUDE
          PATTERN "app" EXCLUDE
          PATTERN "CMakeFiles" EXCLUDE)
