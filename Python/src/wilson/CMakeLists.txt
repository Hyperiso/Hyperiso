add_library(wilson MODULE binding_wilson.cpp)

target_link_libraries(wilson PRIVATE pybind11::module PhysicalModelLib CoreLib)

set_target_properties(wilson PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")