#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_common(py::module &);
void init_core(py::module &);
void init_wilson(py::module &);
void init_observable(py::module &);

PYBIND11_MODULE(pyhyperiso, m) {
    m.doc() = "Python interface for hyperiso project";

    auto common = m.def_submodule("common", "Common functionalities for hyperiso");
    init_common(common);

    auto core = m.def_submodule("core", "Core functionalities for hyperiso");
    init_core(core);

    auto wilson = m.def_submodule("wilson", "Wilson coefficient management");
    init_wilson(wilson);

    auto observable = m.def_submodule("observable", "Observable computation");
    init_observable(observable);
}