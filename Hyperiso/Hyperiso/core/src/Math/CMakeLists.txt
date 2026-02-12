cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Math)

if (NOT TARGET CommonLib)
    add_subdirectory(../Common CommonLib)
endif()

    add_library(MathLib STATIC 
        domain/analysis/contour.cpp
        domain/analysis/diffcalc.cpp
        domain/analysis/functions.cpp
        domain/analysis/gsl_wrappers.cpp
        domain/analysis/integration.cpp
        domain/analysis/optimization.cpp
        domain/linalg/Matrix.cpp
        domain/linalg/scalar.cpp
        domain/linalg/scalar.tpp
        domain/special/exponential_integral.cpp
        domain/special/ObsMath.cpp
        domain/special/polylog.cpp
        domain/special/special_generic.cpp
        domain/special/special_SM.cpp
        domain/special/special_THDM.cpp
        domain/special/special_SUSY.cpp     
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