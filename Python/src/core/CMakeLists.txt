add_library(core MODULE binding_core.cpp)

target_include_directories(core PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(core PRIVATE ${PYTHON_LIBRARIES} pybind11::module CoreLib)

set_target_properties(core PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")