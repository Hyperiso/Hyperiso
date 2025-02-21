#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_utilities(py::module &m) {
}

PYBIND11_MODULE(wilson_interface, m) {
    m.doc() = "Python interface for utility functions";

    init_utilities(m);
}
