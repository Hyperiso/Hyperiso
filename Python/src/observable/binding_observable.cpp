#include <pybind11/pybind11.h>
#include "ObservableInterface.h"

namespace py = pybind11;

void init_observable(py::module &m) {

    // py::enum_<Observables>(m, "Observables", "Defines available observables.")
    //     .value("BR_BS_MUMU", Observables::BR_BS_MUMU, "Branching ratio of Bs -> mu+ mu-")
    //     .value("BR_BS_MUMU_UNTAG", Observables::BR_BS_MUMU_UNTAG, "Untagged branching ratio of Bs -> mu+ mu-")
    //     .value("BR_BD_MUMU", Observables::BR_BD_MUMU, "Branching ratio of Bd -> mu+ mu-")
    //     .value("BR_BU_TAUNU", Observables::BR_BU_TAUNU, "Branching ratio of Bu -> tau nu")
    //     .value("ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA", Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "Isospin asymmetry in B -> K* gamma")
    //     .export_values();
        
    py::class_<ObservableInterface, std::shared_ptr<ObservableInterface>>(m, "ObservableInterface")
        .def(py::init<>())
        .def("compute_observable", &ObservableInterface::compute_observable)
        .def("set_param", &ObservableInterface::set_param);
}