#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "WilsonManager.h"
#include "MartyWilson.h"

namespace py = pybind11;

void init_coefficient_groups(py::module &m) {
    py::class_<CoefficientGroup, std::shared_ptr<CoefficientGroup>>(m, "CoefficientGroup")
        .def("get_matching", &CoefficientGroup::getMatching)
        .def("get_run", &CoefficientGroup::getRun)
        .def("set_q_match", &CoefficientGroup::set_Q_match)
        .def("set_q_run", &CoefficientGroup::set_Q_run);

    py::class_<BCoefficientGroup, CoefficientGroup, std::shared_ptr<BCoefficientGroup>>(m, "BCoefficientGroup")
        .def(py::init<>());

    py::class_<BCoefficientGroupMarty, CoefficientGroup, std::shared_ptr<BCoefficientGroupMarty>>(m, "BCoefficientGroupMarty")
        .def(py::init<>());

    py::class_<BPrimeCoefficientGroup, CoefficientGroup, std::shared_ptr<BPrimeCoefficientGroup>>(m, "BPrimeCoefficientGroup")
        .def(py::init<>());

    py::class_<BPrimeCoefficientGroupMarty, CoefficientGroup, std::shared_ptr<BPrimeCoefficientGroupMarty>>(m, "BPrimeCoefficientGroupMarty")
        .def(py::init<>());

    py::class_<BScalarCoefficientGroup, CoefficientGroup, std::shared_ptr<BScalarCoefficientGroup>>(m, "BScalarCoefficientGroup")
        .def(py::init<>());

    py::class_<BScalarCoefficientGroupMarty, CoefficientGroup, std::shared_ptr<BScalarCoefficientGroupMarty>>(m, "BScalarCoefficientGroupMarty")
        .def(py::init<>());
}

void init_wilson_coefficient(py::module &m) {
    py::class_<WilsonCoefficient, std::shared_ptr<WilsonCoefficient>>(m, "WilsonCoefficient")
        .def("get_coefficient_matching_value", &WilsonCoefficient::get_CoefficientMatchingValue)
        .def("get_coefficient_run_value", &WilsonCoefficient::get_CoefficientRunValue)
        .def("set_coefficient_matching_value", &WilsonCoefficient::set_CoefficientMatchingValue)
        .def("set_coefficient_run_value", &WilsonCoefficient::set_WilsonCoeffRun)
        .def("get_q_match", &WilsonCoefficient::get_Q_match)
        .def("get_q", &WilsonCoefficient::get_Q);
}

void init_coefficient_manager(py::module &m) {
    py::class_<CoefficientManager, std::shared_ptr<CoefficientManager>>(m, "CoefficientManager")
        .def_static("get_instance", [](const std::string &modelName) {
            return CoefficientManager::GetInstance(modelName);
        }, py::return_value_policy::reference)

        .def_static("initialize", &CoefficientManager::initialize, py::arg("lhaFile"), py::arg("models"),
            "Initialise MemoryManager avec un fichier LHA et une liste de modèles.")
        

        .def("register_coefficient_group", &CoefficientManager::registerCoefficientGroup)
        .def("get_state", &CoefficientManager::get_state)
        .def("set_q_match", &CoefficientManager::setQMatch)
        .def("set_group_scale", &CoefficientManager::setGroupScale)
        .def("set_matching_coefficient", &CoefficientManager::setMatchingCoefficient)
        .def("set_run_coefficient", &CoefficientManager::setRunCoefficient)
        .def("get_matching_coefficient", &CoefficientManager::getMatchingCoefficient)
        .def("get_run_coefficient", &CoefficientManager::getRunCoefficient)
        .def("get_coefficient_group", &CoefficientManager::getCoefficientGroup, py::return_value_policy::reference);
}

void init_wilson_parameters(py::module &m) {
    py::class_<Wilson_parameters, std::shared_ptr<Wilson_parameters>>(m, "WilsonParameters")
        .def_static("get_instance", &Wilson_parameters::GetInstance, py::return_value_policy::reference)
        .def("set_mu", &Wilson_parameters::SetMu)
        .def("set_mu_w", &Wilson_parameters::SetMuW)
        .def("set_gen", &Wilson_parameters::set_gen);
}

PYBIND11_MODULE(wilson, m) {
    m.doc() = "Python interface for Wilson coefficient management";

    init_coefficient_groups(m);
    init_wilson_coefficient(m);
    init_coefficient_manager(m);
    init_wilson_parameters(m);
}
