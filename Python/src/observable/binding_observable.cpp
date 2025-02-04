#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "ObservableInterface.h"

namespace py = pybind11;

void init_observable(py::module &m) {
        
    py::class_<ObservableInterface, std::shared_ptr<ObservableInterface>>(m, "ObservableInterface")
        .def(py::init<>())
        .def("add_observable", &ObservableInterface::add_observable)
        .def("add_observables", &ObservableInterface::add_observables)
        .def("add_observable_parameter", &ObservableInterface::add_observable_parameter)
        .def("add_observable_parameters", &ObservableInterface::add_observable_parameters)
        .def("compute_observable", &ObservableInterface::compute_observable)
        .def("compute_all_observables", &ObservableInterface::compute_all_observables)
        .def("compute_uncertainty", &ObservableInterface::compute_uncertainty)
        .def("compute_leading_uncertainties", &ObservableInterface::compute_leading_uncertainties)
        .def("compute_all_uncertainties", &ObservableInterface::compute_all_uncertainties)
        .def("compute_chi2", &ObservableInterface::compute_chi2)
        .def("set_param", &ObservableInterface::set_param)
        .def("get_param", &ObservableInterface::get_param);
}