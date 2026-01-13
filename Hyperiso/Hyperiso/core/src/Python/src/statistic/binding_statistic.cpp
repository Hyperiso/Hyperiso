#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "StatisticInterface.h"

namespace py = pybind11;

void init_statistic(py::module &m) {

    py::class_<StatisticConfig>(m, "StatisticConfig")
    .def(py::init<>())
    .def_readwrite("obss", &StatisticConfig::obss)
    .def_readwrite("p_specs", &StatisticConfig::p_specs)
    .def_readwrite("MC_draws", &StatisticConfig::MC_draws)
    .def_readwrite("skew_abs_threshold", &StatisticConfig::skew_abs_threshold)
    .def_readwrite("MLE_max_iter", &StatisticConfig::MLE_max_iter)
    .def_readwrite("MLE_tol", &StatisticConfig::MLE_tol);

    py::class_<FitResultWithMaps>(m, "FitResultWithMaps");

    py::class_<StatisticInterface, std::shared_ptr<StatisticInterface>>(m, "StatisticInterface")
     .def(py::init<StatisticConfig>(), py::arg("config"))

     .def("compute_uncertainties", &StatisticInterface::compute_uncertainties);
    //  .def("compute_MLE", &StatisticInterface::compute_MLE);
    
}
