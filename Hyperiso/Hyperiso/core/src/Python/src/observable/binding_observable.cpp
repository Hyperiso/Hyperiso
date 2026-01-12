#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "ObservableInterface.h"
// #include "Compound.h"

namespace py = pybind11;

void init_observable(py::module &m) {

     // py::class_<Estimate>(m, "Estimate")
     //      .def_readwrite("central_value", &Estimate::central_value)
     //      .def_readwrite("stat_std", &Estimate::stat_std)
     //      .def_readwrite("syst_std", &Estimate::syst_std)
     //      .def("combined_std", &Estimate::combined_std)
     //      .def("__repr__",
     //           [](const Estimate& e) {
     //                std::ostringstream oss;
     //                oss << e;
     //                return oss.str();
     //           }
     //      );

    py::class_<ObservableInterface, std::shared_ptr<ObservableInterface>>(m, "ObservableInterface")
     .def(py::init<>())

     // add_observable( … )
     .def("add_observable",
          py::overload_cast<Observables, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)
     .def("add_observable",
          py::overload_cast<ObservableId, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)

     // add_observables( … ) – tu exposes seulement certaines surcharges
     .def("add_observables", 
          py::overload_cast<std::map<Observables, QCDOrder>, bool>(&ObservableInterface::add_observables),
          py::arg("obss"), py::arg("add_dependencies") = false)
     .def("add_observables", 
          py::overload_cast<Decays, QCDOrder, bool>(&ObservableInterface::add_observables),
          py::arg("decay"), py::arg("order"), py::arg("add_dependencies") = false)

     .def("add_observable_parameter",
          py::overload_cast<Observables, ParamId>(&ObservableInterface::add_observable_parameter))
     .def("add_observable_parameter",
          py::overload_cast<ObservableId, ParamId>(&ObservableInterface::add_observable_parameter))

     .def("add_observable_parameters",
          py::overload_cast<Observables, std::unordered_set<ParamId>>(&ObservableInterface::add_observable_parameters))
     .def("add_observable_parameters",
          py::overload_cast<ObservableId, std::unordered_set<ParamId>>(&ObservableInterface::add_observable_parameters))

     // compute_observable(…) const
     .def("compute_observable",
          py::overload_cast<Observables>(&ObservableInterface::compute_observable, py::const_),
          py::arg("obs"))
     .def("compute_observable",
          py::overload_cast<ObservableId>(&ObservableInterface::compute_observable, py::const_),
          py::arg("obs"))

     // // compute_uncertainty(…) const
     // .def("compute_uncertainty",
     //      py::overload_cast<Observables, UncertaintyType>(&ObservableInterface::compute_uncertainty, py::const_),
     //      py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)
     // .def("compute_uncertainty",
     //      py::overload_cast<ObservableId, UncertaintyType>(&ObservableInterface::compute_uncertainty, py::const_),
     //      py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)

     // compute_leading_uncertainties(…) const
     // .def("compute_leading_uncertainties",
     //      py::overload_cast<Observables, size_t, UncertaintyType>(
     //           &ObservableInterface::compute_leading_uncertainties, py::const_),
     //      py::arg("obs"), py::arg("n"), py::arg("u_type") = UncertaintyType::COMBINED)
     // .def("compute_leading_uncertainties",
     //      py::overload_cast<ObservableId, size_t, UncertaintyType>(
     //           &ObservableInterface::compute_leading_uncertainties, py::const_),
     //      py::arg("obs"), py::arg("n"), py::arg("u_type") = UncertaintyType::COMBINED)

     // .def("compute_all_uncertainties", &ObservableInterface::compute_all_uncertainties)
     // .def("compute_all", &ObservableInterface::compute_all)
     // .def("compute_chi2", &ObservableInterface::compute_chi2)

     // remove_observable / remove_observables
     .def("remove_observable",
          py::overload_cast<Observables>(&ObservableInterface::remove_observable),
          py::arg("obs"))
     .def("remove_observable",
          py::overload_cast<ObservableId>(&ObservableInterface::remove_observable),
          py::arg("obs"))

     .def("remove_observables",
          py::overload_cast<std::unordered_set<Observables>>(
               &ObservableInterface::remove_observables),
          py::arg("ids"))
     .def("remove_observables",
          py::overload_cast<Decays>(
               &ObservableInterface::remove_observables),
          py::arg("decay"))

     // get_exp_value / get_exp_uncertainty
     .def("get_exp_value",
          py::overload_cast<Observables>(&ObservableInterface::get_exp_value),
          py::arg("obs"))
     .def("get_exp_value",
          py::overload_cast<ObservableId>(&ObservableInterface::get_exp_value),
          py::arg("obs"))

     .def("get_exp_uncertainty",
          py::overload_cast<Observables, UncertaintyType>(&ObservableInterface::get_exp_uncertainty),
          py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)
     .def("get_exp_uncertainty",
          py::overload_cast<ObservableId, UncertaintyType>(&ObservableInterface::get_exp_uncertainty),
          py::arg("obs"), py::arg("u_type") = UncertaintyType::COMBINED)

     .def("get_current_observables", &ObservableInterface::get_current_observables)
     // .def("get_all_exp", &ObservableInterface::get_all_exp)
     .def("set_param", &ObservableInterface::set_param)
     .def("get_param", &ObservableInterface::get_param);

     // // update_gradient
     // .def("update_gradient",
     //      py::overload_cast<Observables>(&ObservableInterface::update_gradient),
     //      py::arg("obs"))
     // .def("update_gradient",
     //      py::overload_cast<ObservableId>(&ObservableInterface::update_gradient),
     //      py::arg("obs"));
}