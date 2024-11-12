#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "MemoryManager.h"
#include "Parameters.h"

namespace py = pybind11;

void init_memory_manager(py::module &m) {
    py::class_<MemoryManager, std::shared_ptr<MemoryManager>>(m, "MemoryManager")
        .def_static("get_instance", &MemoryManager::GetInstance, py::arg("lhaFile") = "DataBase/example.flha", py::arg("models") = std::vector<int>{0}, py::return_value_policy::reference)
        .def("init", &MemoryManager::init)
        .def("get_input_lha_path", &MemoryManager::getInputLhaPath)
        .def("get_data", &MemoryManager::getData);
}

void init_parameters(py::module &m) {
    py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
        .def_static("get_instance", &Parameters::GetInstance, py::arg("modelId") = 0, py::return_value_policy::reference)
        .def("alpha_s", &Parameters::alpha_s)
        .def("running_mass", &Parameters::running_mass)
        .def("set_block_value", &Parameters::setBlockValue)
        .def("get_qcd_masse", &Parameters::get_QCD_masse);
}

void init_model_strategy(py::module &m) {
    py::class_<ModelStrategy, std::shared_ptr<ModelStrategy>>(m, "ModelStrategy");

    py::class_<SMModelStrategy, ModelStrategy, std::shared_ptr<SMModelStrategy>>(m, "SMModelStrategy")
        .def(py::init<>());

    py::class_<SUSYModelStrategy, ModelStrategy, std::shared_ptr<SUSYModelStrategy>>(m, "SUSYModelStrategy")
        .def(py::init<>());

    py::class_<THDMModelStrategy, ModelStrategy, std::shared_ptr<THDMModelStrategy>>(m, "THDMModelStrategy")
        .def(py::init<>());
}

PYBIND11_MODULE(core, m) {
    m.doc() = "Core functionalities for hyperiso";

    init_memory_manager(m);
    init_parameters(m);
    init_model_strategy(m);
}
