#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "General.h"
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
        

    py::enum_<Observables>(m, "Observables")
        .value("BR_BS_MUMU", Observables::BR_BS_MUMU)
        .value("BR_BS_MUMU_UNTAG", Observables::BR_BS_MUMU_UNTAG)
        .value("BR_BD_MUMU", Observables::BR_BD_MUMU)
        .value("R_TAU_NU", Observables::R_TAU_NU)
        .value("BR_BU_TAU_NU", Observables::BR_BU_TAU_NU)
        .value("ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA", Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA)
        .value("BR_B_XS_GAMMA", Observables::BR_B_XS_GAMMA)
        .value("BR_B__D_TAU_NU", Observables::BR_B__D_TAU_NU)
        .value("XI__D_L_NU", Observables::XI__D_L_NU)
        .export_values();

    py::enum_<QCDOrder>(m, "QCDOrder")
        .value("NONE", QCDOrder::NONE)
        .value("LO", QCDOrder::LO)
        .value("NLO", QCDOrder::NLO)
        .value("NNLO", QCDOrder::NNLO)
        .export_values();

    py::enum_<WCoef>(m, "WCoef")
        .value("C1", WCoef::C1)
        .value("C2", WCoef::C2)
        .value("C3", WCoef::C3)
        .value("C4", WCoef::C4)
        .value("C5", WCoef::C5)
        .value("C6", WCoef::C6)
        .value("C7", WCoef::C7)
        .value("C8", WCoef::C8)
        .value("C9", WCoef::C9)
        .value("C10", WCoef::C10)
        .value("CQ1", WCoef::CQ1)
        .value("CQ2", WCoef::CQ2)
        .value("CP1", WCoef::CP1)
        .value("CP2", WCoef::CP2)
        .value("CP3", WCoef::CP3)
        .value("CP4", WCoef::CP4)
        .value("CP5", WCoef::CP5)
        .value("CP6", WCoef::CP6)
        .value("CP7", WCoef::CP7)
        .value("CP8", WCoef::CP8)
        .value("CP9", WCoef::CP9)
        .value("CP10", WCoef::CP10)
        .value("CPQ1", WCoef::CPQ1)
        .value("CPQ2", WCoef::CPQ2)
        .export_values();

    py::enum_<WGroup>(m, "WGroup")
        .value("B", WGroup::B)
        .value("BPrime", WGroup::BPrime)
        .value("BScalar", WGroup::BScalar)
        .export_values();

    py::enum_<BWilsonBasis>(m, "BWilsonBasis")
        .value("STANDARD", BWilsonBasis::STANDARD)
        .value("TRADITIONAL", BWilsonBasis::TRADITIONAL)
        .export_values();

    py::class_<ParamId>(m, "ParamId")
        .def(py::init<ParameterType, std::string, int>())
        .def_readwrite("type", &ParamId::type)
        .def_readwrite("block", &ParamId::block)
        .def_readwrite("code", &ParamId::code);

    py::class_<ObservableMapper, std::shared_ptr<ObservableMapper>>(m, "ObservableMapper")
        .def_static(
            "str",
            &ObservableMapper::str,
            py::arg("obs")
        )
        .def_static(
            "enum_elt",
            &ObservableMapper::enum_elt,
            py::arg("obs")
        )
        .def_static(
            "get_str",
            &ObservableMapper::get_str
        )
        .def_static(
            "get_enum",
            &ObservableMapper::get_enum);

    py::class_<OrderMapper, std::shared_ptr<OrderMapper>>(m, "OrderMapper")
        .def_static(
            "str",
            &OrderMapper::str,
            py::arg("order")
        )
        .def_static(
            "enum_elt",
            &OrderMapper::enum_elt,
            py::arg("order")
        )
        .def_static(
            "get_str",
            &OrderMapper::get_str
        )
        .def_static(
            "get_enum",
            &OrderMapper::get_enum);

    py::class_<GroupMapper, std::shared_ptr<GroupMapper>>(m, "GroupMapper")
        .def_static(
            "str",
            &GroupMapper::str,
            py::arg("group")
        )
        .def_static(
            "enum_elt",
            &GroupMapper::enum_elt,
            py::arg("group")
        )
        .def_static(
            "get_str",
            &GroupMapper::get_str
        )
        .def_static(
            "get_enum",
            &GroupMapper::get_enum);

    py::class_<WCoefMapper, std::shared_ptr<WCoefMapper>>(m, "WCoefMapper")
        .def_static(
            "str",
            &WCoefMapper::str,
            py::arg("coeff")
        )
        .def_static(
            "enum_elt",
            &WCoefMapper::enum_elt,
            py::arg("coeff")
        )
        .def_static(
            "get_str",
            &WCoefMapper::get_str
        )
        .def_static(
            "get_enum",
            &WCoefMapper::get_enum)
        .def_static(
            "get_group",
            &WCoefMapper::get_group,
            py::arg("group")
        );

    py::class_<ParameterTypeMapper, std::shared_ptr<ParameterTypeMapper>>(m, "ParameterTypeMapper")
        .def_static(
            "str",
            &ParameterTypeMapper::str,
            py::arg("type")
        )
        .def_static(
            "enum_elt",
            &ParameterTypeMapper::enum_elt,
            py::arg("type")
        )
        .def_static(
            "get_str",
            &ParameterTypeMapper::get_str
        )
        .def_static(
            "get_enum",
            &ParameterTypeMapper::get_enum);
    
    py::class_<ModelMapper, std::shared_ptr<ModelMapper>>(m, "ModelMapper")
        .def_static(
            "str",
            &ModelMapper::str,
            py::arg("type")
        )
        .def_static(
            "enum_elt",
            &ModelMapper::enum_elt,
            py::arg("type")
        )
        .def_static(
            "get_str",
            &ModelMapper::get_str
        )
        .def_static(
            "get_enum",
            &ModelMapper::get_enum);

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
