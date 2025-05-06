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
    integration.cpp)
set_target_properties(MathLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
find_package(GSL REQUIRED)

set_target_properties(MathLib PROPERTIES PUBLIC_HEADER "Math.h;Math_SUSY.h;Math_THDM.h,Formfactors.h,Matrix.h")

target_include_directories(MathLib PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories(MathLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/MathBSM)

target_link_libraries(MathLib PRIVATE CommonLib ${GSL_LIBRARIES})
