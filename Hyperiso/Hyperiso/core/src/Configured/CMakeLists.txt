# cmake_minimum_required(VERSION 3.20)

# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# project(ConfiguredModule)

if (NOT DEFINED PROJECT_ROOT)
    get_filename_component(PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../" ABSOLUTE)
    if(NOT PROJECT_ROOT MATCHES "/$")
        set(PROJECT_ROOT "${PROJECT_ROOT}/")
    endif()
endif()

if (NOT DEFINED PROJECT_ASSETS_ROOT)
    get_filename_component(PROJECT_ASSETS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../Assets/" ABSOLUTE)
    if(NOT PROJECT_ASSETS_ROOT MATCHES "/$")
        set(PROJECT_ASSETS_ROOT "${PROJECT_ASSETS_ROOT}/")
    endif()
endif()

if (NOT DEFINED PROJECT_TP_ROOT)
    get_filename_component(PROJECT_TP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../Third_party/" ABSOLUTE)
    if(NOT PROJECT_TP_ROOT MATCHES "/$")
        set(PROJECT_TP_ROOT "${PROJECT_TP_ROOT}/")
    endif()
endif()

add_library(ConfiguredModule STATIC dummy.cpp)

configure_file(
    config.hpp.in
    "${CMAKE_BINARY_DIR}/configured_files/include/config.hpp" ESCAPE_QUOTES
)

# target_include_directories(ConfiguredModule INTERFACE "${CMAKE_BINARY_DIR}/configured_files/include")

target_include_directories(ConfiguredModule
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/configured_files/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Configured>
)

include(GNUInstallDirs)

install(TARGETS ConfiguredModule
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Installer les headers générés
install(DIRECTORY ${CMAKE_BINARY_DIR}/configured_files/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Configured
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")