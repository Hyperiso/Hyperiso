cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Common)

add_library(CommonLib STATIC 
    Logger.cpp
    General.cpp
    Observer.cpp
    Mapper/Map.cpp)

set_target_properties(CommonLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
# set_target_properties(CommonLib PROPERTIES PUBLIC_HEADER "Logger.h;General.h;Utils.h")

# target_include_directories(CommonLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_include_directories(CommonLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Mapper)


target_include_directories(CommonLib
  PUBLIC
    # build
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Mapper>
    # install
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Common>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Common/Mapper>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Mapper
)





add_executable(aaaah aaaah.cpp)

target_link_libraries(aaaah CommonLib)

install(TARGETS CommonLib
  EXPORT HyperisoTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Common
  FILES_MATCHING
    PATTERN "*.h"    PATTERN "*.hpp"
    PATTERN "*.tpp"  PATTERN "*.ipp"  PATTERN "*.tcc"  PATTERN "*.inl"
    PATTERN "*.cpp" EXCLUDE
    PATTERN "app"  EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
)