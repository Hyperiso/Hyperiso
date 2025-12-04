#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "StatisticInterface.h"

namespace py = pybind11;

void init_statistic(py::module &m) {

    py::class_<StatisticConfig>(m, "StatisticConfig");
    py::class_<FitResultWithMaps>(m, "FitResultWithMaps");

    py::class_<StatisticInterface, std::shared_ptr<StatisticInterface>>(m, "StatisticInterface")
     .def(py::init<StatisticConfig>(), py::arg("config"))

     .def("compute_uncertainties", &StatisticInterface::compute_uncertainties)
     .def("compute_MLE", &StatisticInterface::compute_MLE);
    
}
