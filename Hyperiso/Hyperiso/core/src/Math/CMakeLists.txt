cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Math)

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()

    add_library(MathLib STATIC 
        scalar.cpp
        bessel_function.cpp
        exponential_integral.cpp
        others_functions.cpp
        polylog.cpp
        wilson_math.cpp
        MathBSM/wilson_susy.cpp
        MathBSM/wilson_THDM.cpp
        linalg/Matrix.cpp
        analysis/integration.cpp
        analysis/root_finding.cpp
        analysis/minimize.cpp
        analysis/differentiate.cpp
        analysis/contour.cpp
        ObsMath.cpp
        )

set_target_properties(MathLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
find_package(GSL REQUIRED)


target_include_directories(MathLib PRIVATE ${GSL_INCLUDE_DIRS})

    target_include_directories(MathLib
    PUBLIC
        # build
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MathBSM>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/linalg>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/analysis>
        # install
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Math>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Math/MathBSM>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Math/linalg>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Math/analysis>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/MathBSM
        ${CMAKE_CURRENT_SOURCE_DIR}/linalg
        ${CMAKE_CURRENT_SOURCE_DIR}/analysis
    )

if (TARGET GSL::gsl)
    target_link_libraries(MathLib PUBLIC CommonLib GSL::gsl)
    else()
    target_include_directories(MathLib PUBLIC ${GSL_INCLUDE_DIRS})
    target_link_libraries(MathLib PUBLIC CommonLib ${GSL_LIBRARIES})
endif()


include(GNUInstallDirs)

install(TARGETS MathLib
    EXPORT HyperisoTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Math
  FILES_MATCHING
    PATTERN "*.h"    PATTERN "*.hpp"
    PATTERN "*.tpp"  PATTERN "*.ipp"  PATTERN "*.tcc"  PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "app"  EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)