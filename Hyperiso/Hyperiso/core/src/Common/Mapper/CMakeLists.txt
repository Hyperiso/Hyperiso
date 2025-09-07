cmake_minimum_required(VERSION 3.16)
project(DynamicRegistryTest LANGUAGES CXX)

# ✅ Adapte ce chemin vers la racine de tes headers/sources existants
#   (là où se trouvent GeneralEnum.h, Map.h/cpp, ObservableMapper.h, etc.)
set(CORE_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/..")  # <-- modifie au besoin

add_executable(test_dynamic_registry
  test_dynamic_registry.cpp
  # On link Map.cpp pour les mapping() si ce n’est pas déjà dans une lib
  ${CORE_SRC_DIR}/Mapper/Map.cpp
  ${CORE_SRC_DIR}/General.cpp
  ${CORE_SRC_DIR}/Logger.cpp
)

target_include_directories(test_dynamic_registry PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
  ${CORE_SRC_DIR}
)

target_compile_features(test_dynamic_registry PRIVATE cxx_std_20)

# Un peu de warnings
if (MSVC)
  target_compile_options(test_dynamic_registry PRIVATE /W4 /permissive-)
else()
  target_compile_options(test_dynamic_registry PRIVATE -Wall -Wextra -Wpedantic)
endif()
