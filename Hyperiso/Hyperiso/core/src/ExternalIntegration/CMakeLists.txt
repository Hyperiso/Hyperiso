cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(ExternalIntegration)

if (NOT DEFINED BUILD_WITH_APP)
    option(BUILD_WITH_APP "Build app/ executable targets" OFF)
endif()

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with bundled SOFTSUSY fallback" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MINUIT2)
    option(BUILD_WITH_MINUIT2 "Build project with Minuit2 support" OFF)
endif()

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()
if (NOT TARGET MathLib)
    add_subdirectory(../Math MathLib)
endif()
if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured Configured)
endif()

add_library(ExternalIntegrationLib STATIC
    Interface.cpp
    2HDMC.cpp
    SoftSusy.cpp
    MinuitCppBackend.cpp
)
set_target_properties(ExternalIntegrationLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(THIRD_PARTY_DIR "" CACHE PATH "Third-party directory (contains MARTY and optional SOFTSUSY)")
if(NOT THIRD_PARTY_DIR)
  get_filename_component(_REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../../.." ABSOLUTE)
  set(THIRD_PARTY_DIR "${_REPO_ROOT}/Third_party")
endif()
message(STATUS "Third-party dir: ${THIRD_PARTY_DIR}")

set(MARTY_ROOT "" CACHE PATH "Existing Marty install prefix (has include/ and lib/)")
set(SOFTSUSY_ROOT "" CACHE PATH "Existing SoftSUSY install prefix or source tree used when BUILD_WITH_SOFTSUSY=ON")

# 2HDMC is vendored and always built as part of ExternalIntegration.
find_package(GSL REQUIRED)
set(HYPERISO_2HDMC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/2HDMC")
set(HYPERISO_2HDMC_SRC "${HYPERISO_2HDMC_ROOT}/src")

if(NOT EXISTS "${HYPERISO_2HDMC_SRC}/THDM.cpp")
    message(FATAL_ERROR "Vendored 2HDMC sources are missing from ${HYPERISO_2HDMC_ROOT}. Re-add Hyperiso/Hyperiso/core/src/ExternalIntegration/third_party/2HDMC.")
endif()

message(STATUS "Building vendored 2HDMC from: ${HYPERISO_2HDMC_ROOT}")
target_sources(ExternalIntegrationLib PRIVATE
    "${HYPERISO_2HDMC_SRC}/THDM.cpp"
    "${HYPERISO_2HDMC_SRC}/SM.cpp"
    "${HYPERISO_2HDMC_SRC}/Util.cpp"
    "${HYPERISO_2HDMC_SRC}/DecayTable.cpp"
    "${HYPERISO_2HDMC_SRC}/Constraints.cpp"
)
target_include_directories(ExternalIntegrationLib PRIVATE
    "${HYPERISO_2HDMC_SRC}"
    ${GSL_INCLUDE_DIRS}
)
if(TARGET GSL::gsl)
    target_link_libraries(ExternalIntegrationLib PRIVATE GSL::gsl)
    if(TARGET GSL::gslcblas)
        target_link_libraries(ExternalIntegrationLib PRIVATE GSL::gslcblas)
    endif()
else()
    target_link_libraries(ExternalIntegrationLib PRIVATE ${GSL_LIBRARIES})
endif()

if (BUILD_WITH_SOFTSUSY)
    add_compile_definitions(BUILD_WITH_SOFTSUSY)
    message(STATUS "Building application with bundled SOFTSUSY fallback")
    include(ExternalProject)

    set(_SOFTSUSY_PREFIX "${THIRD_PARTY_DIR}/SOFTSUSY")
    if(SOFTSUSY_ROOT)
        set(_SOFTSUSY_PREFIX "${SOFTSUSY_ROOT}")
    endif()

    ExternalProject_Add(SOFTSUSY
        GIT_REPOSITORY "https://github.com/BAllanach/softsusy.git"
        GIT_TAG "v4.1.19"
        PREFIX "${_SOFTSUSY_PREFIX}"
        CONFIGURE_COMMAND ""
        COMMAND autoreconf -fiv
        COMMAND ./configure
        BUILD_COMMAND "make"
        BUILD_IN_SOURCE TRUE
        INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(SOFTSUSY SOURCE_DIR)
    target_include_directories(ExternalIntegrationLib PRIVATE "${SOURCE_DIR}")
    add_dependencies(ExternalIntegrationLib SOFTSUSY)
endif()

if (BUILD_WITH_MARTY)
    # add_definitions(-DBUILD_WITH_MARTY)
    add_compile_definitions(BUILD_WITH_MARTY)
    message(STATUS "Building application with MARTY software")

    include(ExternalProject)
  
    # set(MARTY_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../Third_party/MARTY")
    set(MARTY_PATH "${THIRD_PARTY_DIR}/MARTY")
    set(MARTY_INSTALL_PATH "${MARTY_PATH}/install")

    set(MARTY_CORRECT_COMPILER_VERSION "11")  # Remplacez par la version correcte de votre compilateur

    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh" 
        "#!/bin/bash\n"
        "sed -i 's/set(CMAKE_CXX_STANDARD 17)/set(CMAKE_CXX_STANDARD 20)/' ${MARTY_PATH}/src/CMakeLists.txt\n"
        # "sed -i 's/defaultCompilerVersion=7/defaultCompilerVersion=${MARTY_CORRECT_COMPILER_VERSION}/' ${MARTY_PATH}/src/MARTY/.setup_script.sh\n"
        # "sed -i 's/\\([^:]\\)\\bsize_t\\b/\\1std::size_t/g' ${MARTY_PATH}/src/MARTY/csl/src/multipermutation.cpp\n"
        # "sed -i 's/\\([^:]\\)\\bsize_t\\b/\\1std::size_t/g' ${MARTY_PATH}/src/MARTY/csl/include/multipermutation.h\n"
        
    )
    execute_process(COMMAND chmod +x "${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh")

    set(CUSTOM_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh")
    set(INSTALL_PATH "${MARTY_PATH}/MARTY_INSTALL")
    file(WRITE "${MARTY_PATH}/src/run_setup.sh" 
        "#!/bin/bash\n"
        "cd ${MARTY_PATH}\n"
        "script=\"${MARTY_PATH}/MARTY_INSTALL/lib/libmarty.so\"\n"
        "if [ ! -f \"$script\" ]; then\n"
        "echo \" Installation of marty started \"\n"
        "mkdir MARTY_INSTALL\n"
        "cd src/MARTY\n"
        "mkdir build\n"
        "cd build\n"
        "cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${MARTY_PATH}/MARTY_INSTALL\n"
        "sed -i 's/#include \"abstract.h\"/#include <algorithm>\\n#include \"abstract.h\"/g' ${MARTY_PATH}/src/MARTY/src/csl/abreviation.h\n"
        "sed -i 's/#include \"vector.h\"/#include <algorithm>\\n#include \"vector.h\"/g' ${MARTY_PATH}/src/MARTY/src/csl/diagonalization.h\n"
        "sed -i 's/#pragma once/#include <algorithm>\\n#pragma once/g' ${MARTY_PATH}/src/MARTY/src/csl/hardFactor.h\n"
        "sed -i 's/#include <map>/#include <algorithm>\\n#include <map>/g' ${MARTY_PATH}/src/MARTY/src/csl/index.h\n"
        "sed -i 's/#include \"abstract.h\"/#include <type_traits>\\n#include \"parent.h\"\\n#include \"abstract.h\"/g' ${MARTY_PATH}/src/MARTY/src/csl/replace.h\n"
        "sed -i 's/\\([^:]\\)\\bfind_if\\b/\\1std::find_if/g' ${MARTY_PATH}/src/MARTY/src/csl/variableParent.cpp\n"
        "sed -i 's/#include \"abstract.h\"/#include <algorithm>\\n#include \"abstract.h\"/g' ${MARTY_PATH}/src/MARTY/src/csl/variableParent.h\n"
        "sed -i 's/#include <map>/#include <algorithm>\\n#include <map>/g' ${MARTY_PATH}/src/MARTY/src/grafed/core/latexLink.h\n"
        "sed -i 's/#include <map>/#include <algorithm>\\n#include <map>/g' ${MARTY_PATH}/src/MARTY/src/grafed/gui/latexLink.h\n"
        "make\n"
        "make install\n"
        "else\n"
        "echo \"Installation of Marty already done\" \n"
        "fi\n"
        "bashrc_file=\"$HOME/.bashrc\"\n"
        
        "if ! grep -q \"export PATH=\\\$PATH:${INSTALL_PATH}/bin\" \"\$bashrc_file\"; then
        echo \"export PATH=\\\$PATH:${INSTALL_PATH}/bin\" >> \"\$bashrc_file\"
        fi

        if ! grep -q \"export CPATH=\\\$CPATH:${INSTALL_PATH}/include\" \"\$bashrc_file\"; then
            echo \"export CPATH=\\\$CPATH:${INSTALL_PATH}/include\" >> \"\$bashrc_file\"
        fi

        if ! grep -q \"export C_INCLUDE_PATH=\\\$C_INCLUDE_PATH:${INSTALL_PATH}/include\" \"\$bashrc_file\"; then
            echo \"export C_INCLUDE_PATH=\\\$C_INCLUDE_PATH:${INSTALL_PATH}/include\" >> \"\$bashrc_file\"
        fi

        if ! grep -q \"export LIBRARY_PATH=\\\$LIBRARY_PATH:${INSTALL_PATH}/lib\" \"\$bashrc_file\"; then
            echo \"export LIBRARY_PATH=\\\$LIBRARY_PATH:${INSTALL_PATH}/lib\" >> \"\$bashrc_file\"
        fi

        if ! grep -q \"export LD_LIBRARY_PATH=\\\$LD_LIBRARY_PATH:${INSTALL_PATH}/lib\" \"\$bashrc_file\"; then
            echo \"export LD_LIBRARY_PATH=\\\$LD_LIBRARY_PATH:${INSTALL_PATH}/lib\" >> \"\$bashrc_file\"
        fi

        source \"\$bashrc_file\"

        echo \"Environment variable are uptodate in the .bashrc file.\"
    ")

    execute_process(COMMAND chmod +x "${MARTY_PATH}/src/run_setup.sh")

    # set(MARTY_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/marty)    
    ExternalProject_Add(MARTY
        GIT_REPOSITORY "https://github.com/docbrown1955/marty-public.git"
        GIT_TAG 1.6-beta
        PREFIX ${MARTY_PATH}
        # SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/src"
        # BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/build"
        # INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/install"
        CONFIGURE_COMMAND ""
        COMMAND ${MARTY_PATH}/src/run_setup.sh
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        # COMMAND ${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh
        # CMAKE_ARGS
        #     -DCMAKE_INSTALL_PREFIX="${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/install"
        #     -DCMAKE_BUILD_TYPE=Release
        # BUILD_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_UPDATE ON
        LOG_BUILD ON
    )


    ExternalProject_Get_Property(MARTY INSTALL_DIR)
    message(STATUS "${INSTALL_DIR}")
    # target_include_directories(ExternalIntegrationLib PUBLIC "${INSTALL_DIR}/MARTY_INSTALL/include")
    target_link_libraries(ExternalIntegrationLib PUBLIC "${INSTALL_DIR}/MARTY_INSTALL/lib/libmarty.so")
    # target_link_libraries(ExternalIntegrationLib PUBLIC "${INSTALL_DIR}/MARTY_INSTALL/lib/libcsl.so")
    # target_link_libraries(ExternalIntegrationLib PUBLIC "${INSTALL_DIR}/MARTY_INSTALL/lib/libooptools.so")
    add_dependencies(ExternalIntegrationLib MARTY)

endif()
# set_target_properties(ExternalIntegrationLib PROPERTIES PUBLIC_HEADER "2HDMC.h;Interface.h;Softsusy.h")

# target_include_directories(ExternalIntegrationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


target_link_libraries(ExternalIntegrationLib PUBLIC CommonLib)
target_link_libraries(ExternalIntegrationLib PUBLIC ConfiguredModule)
target_link_libraries(ExternalIntegrationLib PUBLIC MathLib)

target_include_directories(ExternalIntegrationLib
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>

    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

include(GNUInstallDirs)

install(TARGETS ExternalIntegrationLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/ExternalIntegration
  FILES_MATCHING
    PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)

# include(FetchContent)

# FetchContent_Declare(minuit_cpp
#   GIT_REPOSITORY https://github.com/doughague/minuit-cpp.git
#   GIT_TAG b5f2558
#   GIT_SHALLOW TRUE
# )


# FetchContent_MakeAvailable(minuit_cpp)

# if(NOT TARGET MinuitCpp)
#   add_library(MinuitCpp INTERFACE)

#   if(TARGET minuit-cpp)
#     target_link_libraries(MinuitCpp INTERFACE minuit-cpp)          # shared
#   elseif(TARGET minuit-cpp-static)
#     target_link_libraries(MinuitCpp INTERFACE minuit-cpp-static)   # static
#   else()
#     message(FATAL_ERROR "Targets minuit-cpp introuvables (minuit-cpp / minuit-cpp-static).")
#   endif()

#   target_include_directories(MinuitCpp INTERFACE
#     ${minuit_cpp_SOURCE_DIR}/include
#   )

#   add_library(MinuitCpp::MinuitCpp ALIAS MinuitCpp)
# endif()

# target_link_libraries(ExternalIntegrationLib PUBLIC MinuitCpp::MinuitCpp)

include(FetchContent)
include(CheckIncludeFile)

FetchContent_Declare(minuit_cpp
  GIT_REPOSITORY https://github.com/doughague/minuit-cpp.git
  GIT_TAG b5f2558
  GIT_SHALLOW TRUE
)

FetchContent_GetProperties(minuit_cpp)
if(NOT minuit_cpp_POPULATED)
  FetchContent_Populate(minuit_cpp)
endif()

check_include_file("getopt.h" HAVE_GETOPT_H)
if(NOT HAVE_GETOPT_H)
  set(HAVE_GETOPT_H 0)
endif()

set(MINUIT_CPP_MONIKER "minuit-cpp")
set(MINUIT_CPP_TITLE "C++ Minimization Package")
set(MINUIT_CPP_AUTHOR "Doug Hague")
set(MINUIT_CPP_YEAR "2013")
set(MINUIT_CPP_VERSION_MAJOR "5")
set(MINUIT_CPP_VERSION_MINOR "29")
set(MINUIT_CPP_VERSION_PATCH "00")
set(MINUIT_CPP_VERSION "${MINUIT_CPP_VERSION_MAJOR}.${MINUIT_CPP_VERSION_MINOR}.${MINUIT_CPP_VERSION_PATCH}")
set(MINUIT_CPP_SOURCE_DIR "${minuit_cpp_SOURCE_DIR}")
set(MINUIT_CPP_BINARY_DIR "${minuit_cpp_BINARY_DIR}")

set(_HYPERISO_SAVED_PROJECT_NAME "${PROJECT_NAME}")
set(PROJECT_NAME "MINUIT_CPP")

file(MAKE_DIRECTORY "${minuit_cpp_BINARY_DIR}/include/minuit-cpp")

configure_file(
  "${minuit_cpp_SOURCE_DIR}/src/Config.hh.in"
  "${minuit_cpp_BINARY_DIR}/include/minuit-cpp/Config.hh"
  @ONLY
)

set(PROJECT_NAME "${_HYPERISO_SAVED_PROJECT_NAME}")

file(GLOB MINUIT_CPP_SOURCES CONFIGURE_DEPENDS
  "${minuit_cpp_SOURCE_DIR}/src/lib/*.cxx"
)

add_library(MinuitCpp STATIC ${MINUIT_CPP_SOURCES})
add_library(MinuitCpp::MinuitCpp ALIAS MinuitCpp)

set_target_properties(MinuitCpp PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)

target_include_directories(MinuitCpp
  PUBLIC
    $<BUILD_INTERFACE:${minuit_cpp_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${minuit_cpp_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(MinuitCpp PUBLIC cxx_std_20)

target_link_libraries(ExternalIntegrationLib PRIVATE MinuitCpp::MinuitCpp)

install(TARGETS MinuitCpp
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)