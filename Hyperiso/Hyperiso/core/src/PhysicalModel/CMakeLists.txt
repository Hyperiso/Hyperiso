cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PhysicalModel)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MINUIT2)
    option(BUILD_WITH_MINUIT2 "Build project with Minuit2 support" OFF)
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
    domain/Helper/MesonMixingRunningParameters.cpp
    domain/Helper/BWilsonRunningParameters.cpp
    domain/Helper/WilsonParametersHelper.cpp
    domain/Wilson.cpp
    domain/MartyWilson.cpp
    domain/WCoef/BWilson.cpp
    domain/WCoef/BPrimeWilson.cpp
    domain/WCoef/BScalarWilson.cpp
    domain/WCoef/KWilson.cpp
    domain/WGroup/BWilsonGroup.cpp
    domain/WGroup/BPrimeWilsonGroup.cpp
    domain/WGroup/BScalarWilsonGroup.cpp
    domain/WGroup/KWilsonGroup.cpp
    domain/WCoef/ChargedCurrentWilson.cpp
    domain/WGroup/ChargedCurrentsWilsonGroup.cpp
    domain/WCoef/DChargedCurrentWilson.cpp
    domain/WGroup/DChargedCurrentsWilsonGroup.cpp
    domain/WCoef/KChargedCurrentWilson.cpp
    domain/WGroup/KChargedCurrentsWilsonGroup.cpp
    domain/WCoef/PIChargedCurrentWilson.cpp
    domain/WGroup/PIChargedCurrentsWilsonGroup.cpp
    domain/WCoef/MesonMixingWilson.cpp
    domain/WGroup/MesonMixingWilsonGroup.cpp
    domain/WilsonGroup.cpp
    domain/WilsonManager.cpp
    domain/MartyWilson.cpp
    domain/GroupDef/GroupDefinitions.cpp
    domain/GroupDef/GroupDef_B.cpp
    domain/GroupDef/GroupDef_CC_bc.cpp
    domain/GroupDef/GroupDef_CC_bu.cpp
    domain/GroupDef/GroupDef_CC_cs.cpp
    domain/GroupDef/GroupDef_CC_cd.cpp
    domain/GroupDef/GroupDef_CC_su.cpp
    domain/GroupDef/GroupDef_CC_du.cpp
    domain/GroupDef/GroupDef_BPrime.cpp
    domain/GroupDef/GroupDef_BScalar.cpp
    domain/GroupDef/GroupDef_MesonMixing.cpp
    domain/GroupDef/GroupDef_K.cpp
    domain/WilsonCoefficientRegistry.cpp
    domain/CoefficientGroupBuilder.cpp
    adapters/ParameterProxy.cpp
    adapters/WilsonParamComposer.cpp
    adapters/ScaleSetter.cpp
    adapters/WilsonBuilder.cpp
    adapters/WilsonProvider.cpp
    adapters/BlockProxy.cpp
    adapters/WilsonInterface.cpp
    domain/Helper/BSM/SUSYParametersHelper.cpp
    domain/Helper/BSM/THDMParametersHelper.cpp
    domain/WCoef/BSM/BWilsonSUSY.cpp
    domain/WCoef/BSM/BWilsonTHDM.cpp
    domain/WCoef/BSM/ChargedCurrentWilsonSUSY.cpp
    domain/WCoef/BSM/ChargedCurrentWilsonTHDM.cpp
    domain/WCoef/BSM/DChargedCurrentWilsonSUSY.cpp
    domain/WCoef/BSM/DChargedCurrentWilsonTHDM.cpp
    domain/WCoef/BSM/KChargedCurrentWilsonSUSY.cpp
    domain/WCoef/BSM/KChargedCurrentWilsonTHDM.cpp
    domain/WCoef/BSM/PIChargedCurrentWilsonSUSY.cpp
    domain/WCoef/BSM/PIChargedCurrentWilsonTHDM.cpp
    )

set_target_properties(PhysicalModelLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
# set_target_properties(PhysicalModelLib PROPERTIES PUBLIC_HEADER "Interpolator.h;WilsonManager.h;Wilson.h;WilsonGroup.h;MartyWilson.h;Wilson_parameters.h;epsilon_calculator.h;susy_parameters.h;thdm_parameters;Wilson_susy.h;Wilson_THDM.h;WilsonManager.h")

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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/GroupDef>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Helper>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/Helper/BSM>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/WCoef>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/WCoef/BSM>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/WGroup>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/domain/WGroup/BSM>

    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/adapters>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/ports>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/GroupDef>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/Helper>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/Helper/BSM>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/WCoef>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/WCoef/BSM>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/WGroup>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/PhysicalModel/domain/WGroup/BSM>
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

add_executable(CUSTOM app/main_custom.cpp)
target_link_libraries(CUSTOM PhysicalModelLib)

add_executable(kaons_test app/main_kaons.cpp)
target_link_libraries(kaons_test PhysicalModelLib)

add_executable(mixing_test app/main_mixing.cpp)
target_link_libraries(mixing_test PhysicalModelLib)

add_executable(charged_test app/main_charged_current.cpp)
target_link_libraries(charged_test PhysicalModelLib)

add_executable(INPUT_WILSON app/main_fwcoef.cpp)
target_link_libraries(INPUT_WILSON PhysicalModelLib)

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
