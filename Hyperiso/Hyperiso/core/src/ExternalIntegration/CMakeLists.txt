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
set_target_properties(ExternalIntegrationLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# set(THIRD_PARTY_DIR "${PROJECT_SOURCE_DIR}/Third_party")

# Permettre de forcer le chemin si besoin :
set(THIRD_PARTY_DIR "" CACHE PATH "Third-party directory (contains MARTY, SOFTSUSY, 2HDMC)")
if(NOT THIRD_PARTY_DIR)
  # Remonter depuis ce dossier jusqu'à la racine du repo (~/hyperiso)
  get_filename_component(_REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../../.." ABSOLUTE)
  set(THIRD_PARTY_DIR "${_REPO_ROOT}/Third_party")
endif()
message(STATUS "Third-party dir: ${THIRD_PARTY_DIR}")

# Si déjà installés ailleurs, on peut pointer dessus pour éviter ExternalProject:
set(MARTY_ROOT "" CACHE PATH "Existing Marty install prefix (has include/ and lib/)")
set(SOFTSUSY_ROOT "" CACHE PATH "Existing SoftSUSY install prefix (has include/ and lib/)")
# (2HDMC: on compile directement les .cpp, mais on peut aussi fournir un src déjà extrait)
set(TWOHDMC_SRC "" CACHE PATH "Existing 2HDMC sources (folder containing src/)")

if (BUILD_WITH_SOFTSUSY)
    # add_definitions(-DBUILD_WITH_SOFTSUSY)    
    add_compile_definitions(BUILD_WITH_SOFTSUSY)
    message(STATUS "Building application with SOFTSUSY software")
    include(ExternalProject)
  
    ExternalProject_Add(SOFTSUSY
        GIT_REPOSITORY "https://github.com/BAllanach/softsusy.git"
        GIT_TAG "v4.1.19"
        PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../Third_party/SOFTSUSY"
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

if (BUILD_WITH_2HDMC)
    # add_definitions(-DBUILD_WITH_2HDMC)
    add_compile_definitions(BUILD_WITH_2HDMC)
    message(STATUS "Building application with 2HDMC software")
    find_package(GSL REQUIRED)

    set(DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/downloads")
    # set(EXTRACT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../Third_party/2HDMC/")
    set(EXTRACT_DIR "${THIRD_PARTY_DIR}/2HDMC")
    set(TAR_FILE "${DOWNLOAD_DIR}/2HDMC-1.8.0.tar.gz")
    set(URL "https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz")

    # Vérifie d'abord si les fichiers sources existent déjà
    if(NOT EXISTS "${EXTRACT_DIR}/src/THDM.cpp")
        message(STATUS "2HDMC sources not found, preparing download and extraction...")

        file(MAKE_DIRECTORY ${DOWNLOAD_DIR})
        file(MAKE_DIRECTORY ${EXTRACT_DIR})

        # Vérifie si le fichier existe mais semble corrompu
        if(EXISTS ${TAR_FILE})
            file(SIZE ${TAR_FILE} TAR_SIZE)
            if(${TAR_SIZE} LESS 10000)
                message(WARNING "2HDMC archive seems corrupted, removing it...")
                file(REMOVE ${TAR_FILE})
            endif()
        endif()

        if(NOT EXISTS ${TAR_FILE})
            message(STATUS "Downloading 2HDMC source archive...")
            execute_process(
                COMMAND wget -O ${TAR_FILE} ${URL}
                WORKING_DIRECTORY ${DOWNLOAD_DIR}
                RESULT_VARIABLE DOWNLOAD_RESULT
            )
            if(NOT DOWNLOAD_RESULT EQUAL 0)
                message(FATAL_ERROR "Failed to download 2HDMC archive from ${URL}")
            endif()
        else()
            message(STATUS "2HDMC archive already exists, skipping download.")
        endif()

        message(STATUS "Extracting 2HDMC source files...")
        execute_process(
            COMMAND tar xzf ${TAR_FILE} --strip-components=1 -C ${EXTRACT_DIR}
            WORKING_DIRECTORY ${DOWNLOAD_DIR}
            RESULT_VARIABLE TAR_RESULT
        )
        if(NOT TAR_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to extract 2HDMC archive.")
        endif()
    else()
        message(STATUS "2HDMC source files already present, skipping download and extraction.")
    endif()

    # Ajout des fichiers à la target
    target_sources(ExternalIntegrationLib PRIVATE
        "${EXTRACT_DIR}/src/THDM.cpp"
        "${EXTRACT_DIR}/src/SM.cpp"
        "${EXTRACT_DIR}/src/Util.cpp"
        "${EXTRACT_DIR}/src/DecayTable.cpp"
    )
    target_link_libraries(ExternalIntegrationLib PRIVATE ${GSL_LIBRARIES})
    target_include_directories(ExternalIntegrationLib PRIVATE "${EXTRACT_DIR}/src" ${GSL_INCLUDE_DIRS})
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

    # Script pour modifier la version du compilateur dans .setup_script.sh
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
set_target_properties(ExternalIntegrationLib PROPERTIES PUBLIC_HEADER "2HDMC.h;Interface.h;Softsusy.h")

# target_include_directories(ExternalIntegrationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


target_link_libraries(ExternalIntegrationLib PUBLIC CommonLib)
target_link_libraries(ExternalIntegrationLib PUBLIC ConfiguredModule)

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