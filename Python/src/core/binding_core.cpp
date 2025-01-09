#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MemoryManager.h"
#include "Parameters.h"

namespace py = pybind11;

void init_core(py::module &m) {

    py::enum_<Model>(m, "Model")
    .value("SM", Model::SM)
    .value("SUSY", Model::SUSY)
    .value("THDM", Model::THDM)
    .value("CUSTOM", Model::CUSTOM)
    .export_values();

    py::enum_<ParameterType>(m, "ParameterType")
        .value("SM", ParameterType::SM)
        .value("SUSY", ParameterType::SUSY)
        .value("THDM", ParameterType::THDM)
        .value("CUSTOM", ParameterType::CUSTOM)
        .value("FLAVOR", ParameterType::FLAVOR)
        .value("WILSON", ParameterType::WILSON)
        .value("FF", ParameterType::FF)
        .export_values();
        
    py::class_<MemoryManager, std::shared_ptr<MemoryManager>>(m, "MemoryManager")
        .def_static("get_instance", &MemoryManager::GetInstance, py::return_value_policy::reference)
        .def("init", &MemoryManager::init,
            py::arg("lhaFile"),
            py::arg_v("model", Model::SM),
            py::arg_v("use_marty", false),
            py::arg("is_spectrum") = false,
            py::arg("has_wilsons") = false,
            py::arg("has_obs") = false)
        .def("get_input_lha_path", &MemoryManager::getInputLhaPath)
        .def("get_data", &MemoryManager::getReader)
        .def("switch_model", &MemoryManager::switch_model)
        .def("switch_lha", &MemoryManager::switch_lha)
        .def("get_blocks_list", &MemoryManager::get_blocks_list)
        .def("get_block_infos", &MemoryManager::get_block_infos);

    py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
        .def_static(
            "get_instance",
            &Parameters::GetInstance,
            py::arg("id") = ParameterType::SM,
            py::return_value_policy::reference
        )
        // .def("alpha_s", &Parameters::alpha_s, py::arg("Q"),
        //      "Compute the strong coupling constant at scale Q.")
        // .def("running_mass", &Parameters::running_mass,
        //      py::arg("quark_mass"), py::arg("q_init"), py::arg("q_end"),
        //      py::arg("option_massb") = "running", py::arg("option_masst") = "pole",
        //      "Compute the running mass of a quark.")
        .def("set_block_value", &Parameters::setBlockValue,
             py::arg("block"), py::arg("code"), py::arg("value"), py::arg("force") = false,
             "Set a value in a specific block.")
        // .def("get_qcd_masse", &Parameters::get_QCD_masse,
        //      py::arg("masstype"),
        //      "Retrieve the QCD mass.")
        .def("exists", &Parameters::exist, py::arg("block"), py::arg("code"),
             "Check if a parameter exists in a block.")
        .def("__call__", [](Parameters &self, const std::string &block, int code) {
                 return self(block, code);
             },
             py::arg("block"), py::arg("code"),
             "Retrieve the value of a parameter in a block.")
        .def_static(
            "get_type",
            &Parameters::GetType,
            py::arg("block"), py::arg("code"),
            "Get the type of a parameter based on its block and code."
        )
        .def("shift_parameter", &Parameters::shiftParameter,
             py::arg("param_id"), py::arg("shift_value"),
             "Shift the value of a parameter by a given value.");
}
