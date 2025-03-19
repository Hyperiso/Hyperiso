cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Core)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()
if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured ConfiguredModule)
endif()
if (NOT TARGET DataBaseLib)
    add_subdirectory(../DataBase DataBaseLib)
endif()
if (NOT TARGET MathLib)
    add_subdirectory(../Math MathLib)
endif()
if (NOT TARGET ExternalIntegrationLib)
    add_subdirectory(../ExternalIntegration ExternalIntegrationLib)
endif()


add_library(CoreLib STATIC 
    adapters/DBMemento.cpp
    adapters/ParamBlockAdapter.cpp
    adapters/APIAdapter.cpp
    adapters/MartyAdapter.cpp
    adapters/HyperisoMaster.cpp
    adapters/CorrelationProvider.cpp
    adapters/ParameterProvider.cpp
    adapters/SpectrumCalculator.cpp
    adapters/CorrelationAdapter.cpp
    adapters/CompositeParamCreator.cpp
    adapters/QCDProvider.cpp
    domain/CorrelationRepo.cpp
    domain/ParameterRouter.cpp
    domain/Parameters.cpp
    domain/Parameter.cpp
    domain/MemoryManager.cpp
    domain/Block.cpp
    domain/BlockAccessor.cpp
    domain/DependentBlockManager.cpp
    infrastructure/QCDHelper.cpp
    )
set_target_properties(CoreLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# set_target_properties(CoreLib PROPERTIES PUBLIC_HEADER "MemoryManager.h;Parameters.h;Parameter.h;QCDHelper.h;BlocksCreator.h;BlockAccessor.h;CorrelationRepo.h;CorrelationCreator.h;ParameterRouter.h;")

target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/adapters)
target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/domain)
target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/infrastructure)
target_include_directories(CoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ports)

target_link_libraries(CoreLib PUBLIC CommonLib DataBaseLib ConfiguredModule ExternalIntegrationLib MathLib)

if (${ENABLE_WARNINGS})
    target_set_warnings(CoreLib
                        ${ENABLE_WARNINGS}
                        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

add_executable(TEST_CORE app/main_core.cpp)
target_link_libraries(TEST_CORE CoreLib)

# add_executable(TEST main_parameters.cpp)
# target_link_libraries(TEST CoreLib)