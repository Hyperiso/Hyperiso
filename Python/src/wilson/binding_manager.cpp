#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "WilsonManager.h"

namespace py = pybind11;

void init_other_manager(py::module &m) {
}

PYBIND11_MODULE(wilson_interface, m) {
    m.doc() = "Python interface for additional managers";

    init_other_manager(m);
}
