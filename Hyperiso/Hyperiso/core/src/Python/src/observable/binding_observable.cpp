#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "ObservableInterface.h"

namespace py = pybind11;

void init_observable(py::module &m) {
    py::class_<ObservableInterface, std::shared_ptr<ObservableInterface>>(m, "ObservableInterface")
        .def(py::init<>()) 
        .def("add_observable", &ObservableInterface::add_observable,
             py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)
        .def("add_observables", 
             py::overload_cast<std::map<Observables, QCDOrder>, bool>(&ObservableInterface::add_observables),
             py::arg("obss"), py::arg("add_dependencies") = false)
        .def("add_observables", 
             py::overload_cast<Decays, QCDOrder, bool>(&ObservableInterface::add_observables),
             py::arg("decay"), py::arg("order"), py::arg("add_dependencies") = false)
        .def("add_observable_parameter", &ObservableInterface::add_observable_parameter)
        .def("add_observable_parameters", &ObservableInterface::add_observable_parameters)
        .def("compute_observable", &ObservableInterface::compute_observable)
        .def("compute_uncertainty", &ObservableInterface::compute_uncertainty,
             py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)
        .def("compute_leading_uncertainties", &ObservableInterface::compute_leading_uncertainties,
             py::arg("obs"), py::arg("n"), py::arg("u_type") = UncertaintyType::COMBINED)
        .def("compute_all_uncertainties", &ObservableInterface::compute_all_uncertainties)
        .def("compute_all", &ObservableInterface::compute_all)
        .def("compute_chi2", &ObservableInterface::compute_chi2)
        .def("remove_observable", &ObservableInterface::remove_observable)
        .def("remove_observables", 
             py::overload_cast<std::unordered_set<Observables>>(&ObservableInterface::remove_observables))
        .def("remove_observables", 
             py::overload_cast<Decays>(&ObservableInterface::remove_observables))
        .def("get_exp_value", &ObservableInterface::get_exp_value)
        .def("get_exp_uncertainty", &ObservableInterface::get_exp_uncertainty,
             py::arg("id"), py::arg("u_type") = UncertaintyType::COMBINED)
        .def("get_current_observables", &ObservableInterface::get_current_observables)
        .def("get_all_exp", &ObservableInterface::get_all_exp)
        .def("set_param", &ObservableInterface::set_param)
        .def("get_param", &ObservableInterface::get_param)
        .def("get_observable_evaluations", &ObservableInterface::get_observable_evaluations)
        .def("update_gradient", &ObservableInterface::update_gradient);
}