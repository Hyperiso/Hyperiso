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

    include(ExternalProject)
  
    # Add SOFTSUSY as an external project
    ExternalProject_Add(SOFTSUSY
        GIT_REPOSITORY "https://github.com/BAllanach/softsusy.git" # Specify the git repository of SOFTSUSY
        GIT_TAG "v4.1.19" # Specify a commit, tag, or branch
        PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/SOFTSUSY"
        # CONFIGURE_COMMAND "autoreconf -fiv && ./configure"
        CONFIGURE_COMMAND ""
        COMMAND autoreconf -fiv
        COMMAND ./configure
        BUILD_COMMAND "make"
        BUILD_IN_SOURCE TRUE
        INSTALL_COMMAND ""
    )

endif()

if (BUILD_WITH_2HDMC)
    add_definitions(-DBUILD_WITH_2HDMC)

    find_package(GSL REQUIRED)

    set(DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/downloads")
    set(EXTRACT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ExternalIntegration/2HDMC/src")
    set(TAR_FILE "${DOWNLOAD_DIR}/2HDMC-1.8.0.tar.gz")
    set(URL "https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz")

    # Téléchargement via execute_process
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
    # Utilisation de la commande shell pour l'extraction
    execute_process(
        COMMAND tar xzf ${TAR_FILE} --strip-components=1 -C ${EXTRACT_DIR}
        WORKING_DIRECTORY ${DOWNLOAD_DIR}
    )

    # file(GLOB THDM "${EXTRACT_DIR}/*.cpp")
    # set(EXCLUDE_FILE "${EXTRACT_DIR}/HBHS.cpp")    

    # if(EXISTS ${EXCLUDE_FILE})
    #     list(REMOVE_ITEM THDM ${EXCLUDE_FILE})
    # else()
    #     message(WARNING "The file ${EXCLUDE_FILE} does not exist.")
    # endif()

    message(STATUS "2HDMC source files: ${THDM}")

    target_sources(ExternalIntegrationLib PRIVATE "${EXTRACT_DIR}/src/THDM.cpp" "${EXTRACT_DIR}/src/SM.cpp" "${EXTRACT_DIR}/src/Util.cpp" "${EXTRACT_DIR}/src/DecayTable.cpp")
    target_include_directories(ExternalIntegrationLib PUBLIC "${EXTRACT_DIR}/src")
endif()



set_target_properties(ExternalIntegrationLib PROPERTIES PUBLIC_HEADER "2HDMC.h;Interface.h;Softsusy.h")

target_include_directories(ExternalIntegrationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


target_link_libraries(ExternalIntegrationLib PUBLIC CommonLib)
target_link_libraries(ExternalIntegrationLib PUBLIC ConfiguredModule)