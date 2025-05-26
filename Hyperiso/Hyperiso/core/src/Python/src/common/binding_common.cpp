#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include "GeneralEnum.h"
#include "EnumMapper.h"
#include "Map.h"
#include "General.h"
#include "Configs.h"

#define BIND_ENUM_MAPPER(cls, type) \
    py::class_<cls, std::shared_ptr<cls>>(m, #cls) \
        .def_static("str", &cls::str, py::arg("value")) \
        .def_static("enum_elt", &cls::enum_elt, py::arg("name")) \
        .def_static("get_str", &cls::get_str) \
        .def_static("get_enum", &cls::get_enum);

namespace py = pybind11;

void init_common(py::module &m) {

    py::enum_<Model>(m, "Model")
    .value("SM", Model::SM)
    .value("SUSY", Model::SUSY)
    .value("THDM", Model::THDM)
    .value("CUSTOM", Model::CUSTOM)
    .export_values();

    py::enum_<ParameterType>(m, "ParameterType")
        .value("SM", ParameterType::SM)
        .value("BSM", ParameterType::BSM)
        .value("FLAVOR", ParameterType::FLAVOR)
        .value("WILSON", ParameterType::WILSON)
        .value("DECAY", ParameterType::DECAY)
        .value("PASSTHROUGH", ParameterType::PASSTHROUGH)
        .value("OBSERVABLE", ParameterType::OBSERVABLE)
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
        .value("A_FB_B__D_TAU_NU", Observables::A_FB_B__D_TAU_NU)
        .value("P_TAU_B__D_TAU_NU", Observables::P_TAU_B__D_TAU_NU)
        .value("R_D", Observables::R_D)
        .value("BR_B__DSTAR_TAU_NU", Observables::BR_B__DSTAR_TAU_NU)
        .value("A_FB_B__DSTAR_TAU_NU", Observables::A_FB_B__DSTAR_TAU_NU)
        .value("P_TAU_B__DSTAR_TAU_NU", Observables::P_TAU_B__DSTAR_TAU_NU)
        .value("P_D_B__DSTAR_TAU_NU", Observables::P_D_B__DSTAR_TAU_NU)
        .value("R_DSTAR", Observables::R_DSTAR)
        .export_values();

    py::enum_<Decays>(m, "Decays")
        .value("B__D_l_nu", Decays::B__D_l_nu)
        .value("B__Dstar_l_nu", Decays::B__Dstar_l_nu)
        .value("B__Kstar", Decays::B__Kstar)
        .value("B__l_l", Decays::B__l_l)
        .value("B__l_nu", Decays::B__l_nu)
        .value("B__Xs", Decays::B__Xs)
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
        .value("CBlnu_A", WCoef::CBlnu_A)
        .value("CBlnu_P", WCoef::CBlnu_P)
        .value("C_V1", WCoef::C_V1)
        .value("C_V2", WCoef::C_V2)
        .value("C_S1", WCoef::C_S1)
        .value("C_S2", WCoef::C_S2)
        .value("C_T", WCoef::C_T)
        .export_values();

    py::enum_<WGroup>(m, "WGroup")
        .value("B", WGroup::B)
        .value("BPrime", WGroup::BPrime)
        .value("BScalar", WGroup::BScalar)
        .value("Blnu", WGroup::Blnu)
        .value("BCLNU", WGroup::BCLNU)
        .export_values();

    py::enum_<WilsonBasis>(m, "WilsonBasis")
        .value("STANDARD", WilsonBasis::B_STANDARD)
        .value("TRADITIONAL", WilsonBasis::B_TRADITIONAL)
        .export_values();

    py::enum_<ContributionType>(m, "ContributionType")
        .value("SM", ContributionType::SM)
        .value("BSM", ContributionType::BSM)
        .value("TOTAL", ContributionType::TOTAL)
        .export_values();


    py::enum_<MassType>(m, "MassType")
        .value("POLE", MassType::POLE)
        .value("MSBAR", MassType::MSBAR)
        .export_values();

    py::enum_<ScaleType>(m, "ScaleType")
        .value("MATCHING", ScaleType::MATCHING)
        .value("HADRONIC", ScaleType::HADRONIC)
        .export_values();


    // py::class_<ParamId>(m, "ParamId")
    //     .def(py::init<ParameterType, std::string, int>())
    //     .def_readwrite("type", &ParamId::type)
    //     .def_readwrite("block", &ParamId::block)
    //     .def_readwrite("code", &ParamId::code);

    

    // py::class_<QCDHelper, std::shared_ptr<QCDHelper>>(m, "QCDHelper")
    //     .def_static(
    //         "mass_b_1S", &QCDHelper::mass_b_1S
    //     );
    
    BIND_ENUM_MAPPER(ObservableMapper, Observables)
    BIND_ENUM_MAPPER(OrderMapper, QCDOrder)
    BIND_ENUM_MAPPER(GroupMapper, WGroup)
    // BIND_ENUM_MAPPER(WCoefMapper, WCoef)
    BIND_ENUM_MAPPER(ParameterTypeMapper, ParameterType)
    BIND_ENUM_MAPPER(ModelMapper, Model)
    BIND_ENUM_MAPPER(WilsonBasisMapper, WilsonBasis)
    BIND_ENUM_MAPPER(ContributionTypeMapper, ContributionType)
    BIND_ENUM_MAPPER(MassTypeMapper, MassType)
    // BIND_ENUM_MAPPER(ScaleTypeMapper, ScaleType)
    // BIND_ENUM_MAPPER(DecayMapper, Decays)

    py::class_<WCoefMapper, std::shared_ptr<WCoefMapper>>(m, "WCoefMapper")
        .def_static("str", &WCoefMapper::str, py::arg("coef"))
        .def_static("enum_elt", &WCoefMapper::enum_elt, py::arg("name"))
        .def_static("get_str", &WCoefMapper::get_str)
        .def_static("get_enum", &WCoefMapper::get_enum)
        .def_static("get_group", &WCoefMapper::get_group, py::arg("group"))
        .def_static("flha_base", &WCoefMapper::flha_base, py::arg("coef"))
        .def_static("flha_full", &WCoefMapper::flha_full, py::arg("coef"), py::arg("order"), py::arg("type"))
        .def_static("from_flha", &WCoefMapper::from_flha, py::arg("content"), py::arg("structure"));

    py::class_<ScaleTypeMapper, std::shared_ptr<ScaleTypeMapper>>(m, "ScaleTypeMapper")
        .def_static("str", &ScaleTypeMapper::str, py::arg("type"))
        .def_static("enum_elt", &ScaleTypeMapper::enum_elt, py::arg("name"))
        .def_static("get_str", &ScaleTypeMapper::get_str)
        .def_static("get_enum", &ScaleTypeMapper::get_enum)
        .def_static("block", &ScaleTypeMapper::block, py::arg("type"));

    py::class_<DecayMapper, std::shared_ptr<DecayMapper>>(m, "DecayMapper")
        .def_static("str", &DecayMapper::str, py::arg("type"))
        .def_static("enum_elt", &DecayMapper::enum_elt, py::arg("name"))
        .def_static("get_str", &DecayMapper::get_str)
        .def_static("get_enum", &DecayMapper::get_enum)
        .def_static("get_observables", &DecayMapper::get_observables, py::arg("decay"))
        .def_static("get_decay", &DecayMapper::get_decay, py::arg("observable"));

 
    py::class_<LhaID>(m, "LhaID")
        .def(py::init<const std::string&>())
        .def(py::init<long>())
        .def(py::init<const std::vector<long>&>())
        .def("__int__", [](const LhaID& self) { return static_cast<long>(self); })
        .def("to_string", &LhaID::to_string)
        .def("get_parts", &LhaID::get_parts);

    py::class_<ParamId>(m, "ParamId")
        .def(py::init<>())
        .def(py::init<const std::string&, const LhaID&>())
        .def(py::init<ParameterType, const std::string&, const LhaID&>())
        .def("set_parameter_type", &ParamId::set_parameter_type)
        .def_readwrite("type", &ParamId::type)
        .def_readwrite("block", &ParamId::block)
        .def_readwrite("code", &ParamId::code);

    py::class_<DependenciesHelper, std::shared_ptr<DependenciesHelper>>(m, "DependenciesHelper")
        .def_static("get_allowed_parameters", &DependenciesHelper::get_allowed_parameters, py::arg("obs"))
        .def_static("is_param_allowed", &DependenciesHelper::is_param_allowed, py::arg("obs"), py::arg("param"));

    py::class_<LhaParamsHelper, std::shared_ptr<LhaParamsHelper>>(m, "LhaParamsHelper")
        .def_static("get_minimal_content", &LhaParamsHelper::get_minimal_content, py::arg("block_name"));


    py::class_<WilsonBuildConfig>(m, "WilsonBuildConfig")
        .def(py::init<>())
        .def_readwrite("groups", &WilsonBuildConfig::groups)
        .def_readwrite("matching_scale", &WilsonBuildConfig::matching_scale)
        .def_readwrite("hadronic_scale", &WilsonBuildConfig::hadronic_scale)
        .def_readwrite("order", &WilsonBuildConfig::order);
    
    py::class_<WilsonRequest>(m, "WilsonRequest")
        .def(py::init<>())
        .def_readwrite("group", &WilsonRequest::group)
        .def_readwrite("coefficient", &WilsonRequest::coefficient)
        .def_readwrite("order", &WilsonRequest::order)
        .def_readwrite("contribution", &WilsonRequest::contribution)
        .def_readwrite("scale_type", &WilsonRequest::scale_type)
        .def_readwrite("sum_qcd_orders", &WilsonRequest::sum_qcd_orders);
    
    py::class_<AlphasConfig>(m, "AlphasConfig")
        .def(py::init<double, MassType, MassType>(), py::arg("scale"), py::arg("m_b_type"), py::arg("m_t_type"))
        .def_readwrite("scale", &AlphasConfig::scale)
        .def_readwrite("m_b_type", &AlphasConfig::m_b_type)
        .def_readwrite("m_t_type", &AlphasConfig::m_t_type);
    
    py::class_<MassConfig, AlphasConfig>(m, "MassConfig")
        .def(py::init<int, double, MassType, MassType>(),
             py::arg("pdg_id"), py::arg("scale"), py::arg("m_b_type"), py::arg("m_t_type"))
        .def_readwrite("pdg_id", &MassConfig::pdg_id);

}
