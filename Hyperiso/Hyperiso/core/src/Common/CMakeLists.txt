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
set_target_properties(CommonLib PROPERTIES PUBLIC_HEADER "Logger.h;General.h;Utils.h")

# target_include_directories(CommonLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_include_directories(CommonLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Mapper)

# target_include_directories(CommonLib
#   PUBLIC
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>               # pendant le build
#     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>           # après install
#   PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}
#     ${CMAKE_CURRENT_SOURCE_DIR}/Mapper
# )

target_include_directories(CommonLib
  PUBLIC
    # Pendant le build, expose directement les dossiers qui contiennent les .h publics
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Mapper>
    # Après install, on les retrouvera là :
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Common>
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

# On copie les headers publics depuis src/Common -> include/hyperiso/Common
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hyperiso/Common
        FILES_MATCHING
          PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inl"
          PATTERN "*.cpp" EXCLUDE
          PATTERN "app" EXCLUDE
          PATTERN "CMakeFiles" EXCLUDE)