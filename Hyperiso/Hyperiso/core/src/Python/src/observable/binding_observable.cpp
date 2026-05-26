#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

#include "ObservableInterface.h"

namespace py = pybind11;

void init_observable(py::module &m) {

    
    py::class_<ObservableValue>(m, "ObservableValue")
        .def(py::init<ObservableId, double>(),
             py::arg("id"), py::arg("value"))
        .def(py::init<ObservableId, double, std::pair<double,double>>(),
             py::arg("id"), py::arg("value"), py::arg("bin"))

        .def_readwrite("id", &ObservableValue::id)
        .def_readwrite("value", &ObservableValue::value)

        .def_property(
            "bin",
            [](const ObservableValue &ov) -> py::object {
                if (ov.bin) return py::cast(*ov.bin);
                return py::none();
            },
            [](ObservableValue &ov, py::object b) {
                if (b.is_none()) ov.bin.reset();
                else ov.bin = b.cast<std::pair<double,double>>();
            })

        .def("__repr__", [](const ObservableValue &ov) {
            std::ostringstream oss;
            oss << "ObservableValue(id=";
            oss << py::str(py::cast(ov.id)).cast<std::string>();
            oss << ", value=" << ov.value;
            if (ov.bin) oss << ", bin=(" << ov.bin->first << ", " << ov.bin->second << ")";
            oss << ")";
            return oss.str();
        });

    py::class_<ObservableInterface, std::shared_ptr<ObservableInterface>>(m, "ObservableInterface")
     .def(py::init<>())

     // add_observable( … )
     .def("add_observable",
          py::overload_cast<Observables, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)
     .def("add_observable",
          py::overload_cast<ObservableId, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)
     .def("add_observable",
          py::overload_cast<BinnedObservableId, QCDOrder, bool>(&ObservableInterface::add_observable),
          py::arg("obs"), py::arg("order"), py::arg("add_dependencies") = false)

     // add_observables( … ) – exposes all public overloads except custom-decay/std::any APIs
     .def("add_observables", 
          py::overload_cast<std::map<Observables, QCDOrder>, bool>(&ObservableInterface::add_observables),
          py::arg("obss"), py::arg("add_dependencies") = false)
     .def("add_observables", 
          py::overload_cast<std::map<ObservableId, QCDOrder>, bool>(&ObservableInterface::add_observables),
          py::arg("obss"), py::arg("add_dependencies") = false)
     .def("add_observables", 
          py::overload_cast<Decays, QCDOrder, bool>(&ObservableInterface::add_observables),
          py::arg("decay"), py::arg("order"), py::arg("add_dependencies") = false)

     .def("add_observable_parameter",
          py::overload_cast<Observables, ParamId>(&ObservableInterface::add_observable_parameter),
          py::arg("obs"), py::arg("pid"))
     .def("add_observable_parameter",
          py::overload_cast<ObservableId, ParamId>(&ObservableInterface::add_observable_parameter),
          py::arg("obs"), py::arg("pid"))

     .def("add_observable_parameters",
          py::overload_cast<Observables, std::unordered_set<ParamId>>(&ObservableInterface::add_observable_parameters),
          py::arg("obs"), py::arg("pids"))
     .def("add_observable_parameters",
          py::overload_cast<ObservableId, std::unordered_set<ParamId>>(&ObservableInterface::add_observable_parameters),
          py::arg("obs"), py::arg("pids"))

     // compute_observable(…) const
     .def("compute_observable",
          py::overload_cast<Observables>(&ObservableInterface::compute_observable, py::const_),
          py::arg("obs"))
     .def("compute_observable",
          py::overload_cast<ObservableId>(&ObservableInterface::compute_observable, py::const_),
          py::arg("obs"))
     .def("compute_all", &ObservableInterface::compute_all)

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
     .def("get_all_ops_deps",
          py::overload_cast<ObservableId>(&ObservableInterface::get_all_ops_deps),
          py::arg("obs"))
     .def("get_all_ops_deps",
          py::overload_cast<Observables>(&ObservableInterface::get_all_ops_deps),
          py::arg("obs"))
     .def("set_param", &ObservableInterface::set_param,
          py::arg("block"), py::arg("code"), py::arg("value"), py::arg("type"))
     .def("get_param", &ObservableInterface::get_param,
          py::arg("block"), py::arg("code"), py::arg("type"))
     .def("reload_params", &ObservableInterface::reload_params)
     .def("enable_obs", &ObservableInterface::enable_obs)
     .def("set_bkstarll_threads", &ObservableInterface::set_bkstarll_threads,
          py::arg("n_threads"))
     .def("set_bkll_threads", &ObservableInterface::set_bkll_threads,
          py::arg("n_threads"))
     .def("set_bsphi_threads", &ObservableInterface::set_bsphi_threads,
          py::arg("n_threads"));

}