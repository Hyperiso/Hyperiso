cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(ExternalIntegration)

if (NOT DEFINED BUILD_WITH_SOFTSUSY)
    option(BUILD_WITH_SOFTSUSY "Build project with SOFTSUSY support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_2HDMC)
    option(BUILD_WITH_2HDMC "Build project with 2HDMC support" OFF)
endif()
if (NOT DEFINED BUILD_WITH_MARTY)
    option(BUILD_WITH_MARTY "Build project with MARTY support" OFF)
endif()

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()

if (NOT TARGET ConfiguredModule)
    add_subdirectory(../Configured Configured)
endif()

add_library(ExternalIntegrationLib STATIC 
    Interface.cpp
    2HDMC.cpp
    SoftSusy.cpp)

if (BUILD_WITH_SOFTSUSY)
    add_definitions(-DBUILD_WITH_SOFTSUSY)
    message(STATUS "Building application with SOFTSUSY software")
    include(ExternalProject)
  
    ExternalProject_Add(SOFTSUSY
        GIT_REPOSITORY "https://github.com/BAllanach/softsusy.git"
        GIT_TAG "v4.1.19"
        PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/SOFTSUSY"
        CONFIGURE_COMMAND ""
        COMMAND autoreconf -fiv
        COMMAND ./configure
        BUILD_COMMAND "make"
        BUILD_IN_SOURCE TRUE
        INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(SOFTSUSY SOURCE_DIR)
    target_include_directories(ExternalIntegrationLib PUBLIC "${SOURCE_DIR}")
    add_dependencies(ExternalIntegrationLib SOFTSUSY)

endif()

if (BUILD_WITH_2HDMC)
    add_definitions(-DBUILD_WITH_2HDMC)
    message(STATUS "Building application with 2HDMC software")
    find_package(GSL REQUIRED)

    set(DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/downloads")
    set(EXTRACT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/2HDMC/src")
    set(TAR_FILE "${DOWNLOAD_DIR}/2HDMC-1.8.0.tar.gz")
    set(URL "https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz")

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DOWNLOAD_DIR}
    )
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E make_directory ${EXTRACT_DIR}
    )
    execute_process(
        COMMAND wget -O ${TAR_FILE} ${URL}
        WORKING_DIRECTORY ${DOWNLOAD_DIR}
    )
    execute_process(
        COMMAND tar xzf ${TAR_FILE} --strip-components=1 -C ${EXTRACT_DIR}
        WORKING_DIRECTORY ${DOWNLOAD_DIR}
    )

    message(STATUS "2HDMC source files: ${THDM}")

    target_sources(ExternalIntegrationLib PRIVATE "${EXTRACT_DIR}/src/THDM.cpp" "${EXTRACT_DIR}/src/SM.cpp" "${EXTRACT_DIR}/src/Util.cpp" "${EXTRACT_DIR}/src/DecayTable.cpp")
    target_link_libraries(ExternalIntegrationLib PRIVATE ${GSL_LIBRARIES})
    target_include_directories(ExternalIntegrationLib PRIVATE "${EXTRACT_DIR}/src" ${GSL_INCLUDE_DIRS})
    
endif()

if (BUILD_WITH_MARTY)
    add_definitions(-DBUILD_WITH_MARTY)
    message(STATUS "Building application with MARTY software")

    include(ExternalProject)
  
    set(MARTY_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY")
    set(MARTY_INSTALL_PATH "${MARTY_PATH}/install")

    set(MARTY_CORRECT_COMPILER_VERSION "11")  # Remplacez par la version correcte de votre compilateur

    # Script pour modifier la version du compilateur dans .setup_script.sh
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh" 
        "#!/bin/bash\n"
        "sed -i 's/defaultCompilerVersion=7/defaultCompilerVersion=${MARTY_CORRECT_COMPILER_VERSION}/' ${MARTY_PATH}/src/MARTY/.setup_script.sh\n"
        "sed -i 's/\\([^:]\\)\\bsize_t\\b/\\1std::size_t/g' ${MARTY_PATH}/src/MARTY/csl/src/multipermutation.cpp\n"
        "sed -i 's/\\([^:]\\)\\bsize_t\\b/\\1std::size_t/g' ${MARTY_PATH}/src/MARTY/csl/include/multipermutation.h\n"
        
    )
    execute_process(COMMAND chmod +x "${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh")


    # Script pour exécuter setup.sh
    file(WRITE "${MARTY_PATH}/src/run_setup.sh" 
        "#!/bin/bash\n"
        "cd ${MARTY_PATH}/src/MARTY\n"
        "source ./setup.sh ${MARTY_INSTALL_PATH}\n"
        "make\n"
        "make install"
    )
    execute_process(COMMAND chmod +x "${MARTY_PATH}/src/run_setup.sh")

        
    ExternalProject_Add(MARTY
        GIT_REPOSITORY "https://github.com/docbrown1955/marty-public.git"
        GIT_TAG "v1.5"
        PREFIX "${MARTY_PATH}"
        # SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/src"
        # BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/MARTY/build"
        CONFIGURE_COMMAND ""
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/MARTY/modify_compiler_version.sh
        COMMAND ${MARTY_PATH}/src/run_setup.sh
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_IN_SOURCE=TRUE
        LOG_DOWNLOAD ON
        LOG_UPDATE ON
        LOG_BUILD ON
    )

    ExternalProject_Get_Property(MARTY INSTALL_DIR)
    target_include_directories(ExternalIntegrationLib PUBLIC "${INSTALL_DIR}/include")
    target_link_libraries(ExternalIntegrationLib PUBLIC "${INSTALL_DIR}/lib/libmarty.so")
    add_dependencies(ExternalIntegrationLib MARTY)

endif()
set_target_properties(ExternalIntegrationLib PROPERTIES PUBLIC_HEADER "2HDMC.h;Interface.h;Softsusy.h")

target_include_directories(ExternalIntegrationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


target_link_libraries(ExternalIntegrationLib PUBLIC CommonLib)
target_link_libraries(ExternalIntegrationLib PUBLIC ConfiguredModule)