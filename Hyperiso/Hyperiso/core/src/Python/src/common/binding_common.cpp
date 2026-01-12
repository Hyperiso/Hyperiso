#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include "GeneralEnum.h"
#include "EnumMapper.h"
#include "Map.h"
#include "Include.h"
#include "Configs.h"

using DH = DependenciesHelper;

// #define BIND_ENUM_MAPPER(cls, type) \
//     py::class_<cls, std::shared_ptr<cls>>(m, #cls) \
//         .def_static("str", &cls::str, py::arg("value")) \
//         .def_static("enum_elt", &cls::enum_elt, py::arg("name")) \
//         .def_static("get_str", &cls::get_str) \
//         .def_static("get_enum", &cls::get_enum);

#define BIND_ENUM_MAPPER(cls, EnumT)                                          \
    py::class_<cls, std::shared_ptr<cls>>(m, #cls)                            \
        /* str(EnumT) -> canonical string */                                  \
        .def_static("str",                                                    \
            py::overload_cast<EnumT>(&cls::str), py::arg("value"))            \
        /* enum_elt(name) -> EnumT (legacy-like) */                           \
        .def_static("enum_elt", &cls::enum_elt_legacy, py::arg("name"))       \
        /* builtins lists (comme avant) */                                    \
        .def_static("get_str",  &cls::get_str)                                \
        .def_static("get_enum", &cls::get_enum)                               \
        /* runtime: builtins + customs */                                     \
        .def_static("get_str_all", &cls::get_str_all);

namespace py = pybind11;

void init_common(py::module &m) {

    py::enum_<Model>(m, "Model")
    .value("SM", Model::SM)
    .value("SUSY", Model::SUSY)
    .value("THDM", Model::THDM)
    .value("MARTY", Model::MARTY)
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
        .value("BR_B__KSTAR_GAMMA", Observables::BR_B__KSTAR_GAMMA)
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
        .value("BR_B__Xs_mu_mu__LOW_Q2", Observables::BR_B__Xs_mu_mu__LOW_Q2)
        .value("BR_B__Xs_mu_mu__HIGH_Q2", Observables::BR_B__Xs_mu_mu__HIGH_Q2)
        .value("BR_B__Xs_tau_tau__HIGH_Q2", Observables::BR_B__Xs_tau_tau__HIGH_Q2)

        .value("DGAMMA_DQ2_B__KSTAR_L_L", Observables::DGAMMA_DQ2_B__KSTAR_L_L)
        .value("DGAMMA_BAR_DQ2_B__KSTAR_L_L", Observables::DGAMMA_BAR_DQ2_B__KSTAR_L_L)
        .value("A_FB_B__KSTAR_L_L", Observables::A_FB_B__KSTAR_L_L)
        .value("Q0_A_FB_B__KSTAR_L_L", Observables::Q0_A_FB_B__KSTAR_L_L)
        .value("A_CP_B__KSTAR_L_L", Observables::A_CP_B__KSTAR_L_L)
        .value("F_L_B__KSTAR_L_L", Observables::F_L_B__KSTAR_L_L)
        .value("F_T_B__KSTAR_L_L", Observables::F_T_B__KSTAR_L_L)
        .value("A_T_1_B__KSTAR_L_L", Observables::A_T_1_B__KSTAR_L_L)
        .value("A_T_2_B__KSTAR_L_L", Observables::A_T_2_B__KSTAR_L_L)
        .value("A_T_3_B__KSTAR_L_L", Observables::A_T_3_B__KSTAR_L_L)
        .value("A_T_4_B__KSTAR_L_L", Observables::A_T_4_B__KSTAR_L_L)
        .value("A_T_5_B__KSTAR_L_L", Observables::A_T_5_B__KSTAR_L_L)
        .value("A_T_RE_B__KSTAR_L_L", Observables::A_T_RE_B__KSTAR_L_L)
        .value("A_T_RE_CPV_B__KSTAR_L_L", Observables::A_T_RE_CPV_B__KSTAR_L_L)
        .value("A_IM_B__KSTAR_L_L", Observables::A_IM_B__KSTAR_L_L)
        .value("ALPHA_K_B__KSTAR_L_L", Observables::ALPHA_K_B__KSTAR_L_L)
        .value("H_T_1_B__KSTAR_L_L", Observables::H_T_1_B__KSTAR_L_L)
        .value("H_T_2_B__KSTAR_L_L", Observables::H_T_2_B__KSTAR_L_L)
        .value("H_T_3_B__KSTAR_L_L", Observables::H_T_3_B__KSTAR_L_L)
        .value("P_2_B__KSTAR_L_L", Observables::P_2_B__KSTAR_L_L)
        .value("P_3_B__KSTAR_L_L", Observables::P_3_B__KSTAR_L_L)
        .value("P_6_B__KSTAR_L_L", Observables::P_6_B__KSTAR_L_L)
        .value("P_8_B__KSTAR_L_L", Observables::P_8_B__KSTAR_L_L)
        .value("P_PRIME_4_B__KSTAR_L_L", Observables::P_PRIME_4_B__KSTAR_L_L)
        .value("P_PRIME_5_B__KSTAR_L_L", Observables::P_PRIME_5_B__KSTAR_L_L)
        .value("P_PRIME_6_B__KSTAR_L_L", Observables::P_PRIME_6_B__KSTAR_L_L)
        .value("P_PRIME_8_B__KSTAR_L_L", Observables::P_PRIME_8_B__KSTAR_L_L)
        .value("S_3_B__KSTAR_L_L", Observables::S_3_B__KSTAR_L_L)
        .value("S_4_B__KSTAR_L_L", Observables::S_4_B__KSTAR_L_L)
        .value("S_5_B__KSTAR_L_L", Observables::S_5_B__KSTAR_L_L)
        .value("S_6C_B__KSTAR_L_L", Observables::S_6C_B__KSTAR_L_L)
        .value("S_7_B__KSTAR_L_L", Observables::S_7_B__KSTAR_L_L)
        .value("S_8_B__KSTAR_L_L", Observables::S_8_B__KSTAR_L_L)
        .value("S_9_B__KSTAR_L_L", Observables::S_9_B__KSTAR_L_L)
        .value("A_3_B__KSTAR_L_L", Observables::A_3_B__KSTAR_L_L)
        .value("A_4_B__KSTAR_L_L", Observables::A_4_B__KSTAR_L_L)
        .value("A_5_B__KSTAR_L_L", Observables::A_5_B__KSTAR_L_L)
        .value("A_6S_B__KSTAR_L_L", Observables::A_6S_B__KSTAR_L_L)
        .value("A_7_B__KSTAR_L_L", Observables::A_7_B__KSTAR_L_L)
        .value("A_8_B__KSTAR_L_L", Observables::A_8_B__KSTAR_L_L)
        .value("A_9_B__KSTAR_L_L", Observables::A_9_B__KSTAR_L_L)
        .value("P_1_CPV_B__KSTAR_L_L", Observables::P_1_CPV_B__KSTAR_L_L)
        .value("P_2_CPV_B__KSTAR_L_L", Observables::P_2_CPV_B__KSTAR_L_L)
        .value("P_3_CPV_B__KSTAR_L_L", Observables::P_3_CPV_B__KSTAR_L_L)
        .value("P_PRIME_4_CPV_B__KSTAR_L_L", Observables::P_PRIME_4_CPV_B__KSTAR_L_L)
        .value("P_PRIME_5_CPV_B__KSTAR_L_L", Observables::P_PRIME_5_CPV_B__KSTAR_L_L)
        .value("P_PRIME_6_CPV_B__KSTAR_L_L", Observables::P_PRIME_6_CPV_B__KSTAR_L_L)
        .value("P_PRIME_8_CPV_B__KSTAR_L_L", Observables::P_PRIME_8_CPV_B__KSTAR_L_L)

        .value("DGAMMA_DQ2_BS__PHI_L_L", Observables::DGAMMA_DQ2_BS__PHI_L_L)
        .value("DGAMMA_BAR_DQ2_BS__PHI_L_L", Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L)
        .value("A_FB_CPV_BS__PHI_L_L", Observables::A_FB_CPV_BS__PHI_L_L)
        .value("F_L_BS_PHI_L_L", Observables::F_L_BS_PHI_L_L)
        .value("A_T_2_BS_PHI_L_L", Observables::A_T_2_BS_PHI_L_L)
        .value("A_T_RE_CPV_BS_PHI_L_L", Observables::A_T_RE_CPV_BS_PHI_L_L)
        .value("A_T_IM_CPV_BS_PHI_L_L", Observables::A_T_IM_CPV_BS_PHI_L_L)
        .value("P_PRIME_4_BS_PHI_L_L", Observables::P_PRIME_4_BS_PHI_L_L)
        .value("P_PRIME_6_BS_PHI_L_L", Observables::P_PRIME_6_BS_PHI_L_L)
        .value("S_2S_BS_PHI_L_L", Observables::S_2S_BS_PHI_L_L)
        .value("S_3_BS_PHI_L_L", Observables::S_3_BS_PHI_L_L)
        .value("S_4_BS_PHI_L_L", Observables::S_4_BS_PHI_L_L)
        .value("S_7_BS_PHI_L_L", Observables::S_7_BS_PHI_L_L)
        .value("A_5_BS_PHI_L_L", Observables::A_5_BS_PHI_L_L)
        .value("A_6C_BS_PHI_L_L", Observables::A_6C_BS_PHI_L_L)
        .value("A_8_BS_PHI_L_L", Observables::A_8_BS_PHI_L_L)
        .value("A_9_BS_PHI_L_L", Observables::A_9_BS_PHI_L_L)
        .value("P_2_CPV_BS_PHI_L_L", Observables::P_2_CPV_BS_PHI_L_L)
        .value("P_3_CPV_BS_PHI_L_L", Observables::P_3_CPV_BS_PHI_L_L)
        .value("P_PRIME_5_CPV_BS_PHI_L_L", Observables::P_PRIME_5_CPV_BS_PHI_L_L)
        .value("P_PRIME_8_CPV_BS_PHI_L_L", Observables::P_PRIME_8_CPV_BS_PHI_L_L)
        .value("Q_8M_BS_PHI_L_L", Observables::Q_8M_BS_PHI_L_L)
        .value("Q_8P_BS_PHI_L_L", Observables::Q_8P_BS_PHI_L_L)
        .value("Q_9_BS_PHI_L_L", Observables::Q_9_BS_PHI_L_L)

        .value("DGAMMA_DQ2_B__K_L_L", Observables::DGAMMA_DQ2_B__K_L_L)
        .value("A_FB_B__K_L_L", Observables::A_FB_B__K_L_L)
        .value("F_H_B__K_L_L", Observables::F_H_B__K_L_L)

        .value("DGAMMA_DQ2_CP_AVG_LAMBDA_B__LAMBDA_L_L", Observables::DGAMMA_DQ2_CP_AVG_LAMBDA_B__LAMBDA_L_L)
        .value("A_FB_L_LAMBDA_B__LAMBDA_L_L", Observables::A_FB_L_LAMBDA_B__LAMBDA_L_L)
        .value("A_FB_H_LAMBDA_B__LAMBDA_L_L", Observables::A_FB_H_LAMBDA_B__LAMBDA_L_L)
        .value("A_FB_LH_LAMBDA_B__LAMBDA_L_L", Observables::A_FB_LH_LAMBDA_B__LAMBDA_L_L)
        .value("F_L_LAMBDA_B__LAMBDA_L_L", Observables::F_L_LAMBDA_B__LAMBDA_L_L)
        .value("F_T_LAMBDA_B__LAMBDA_L_L", Observables::F_T_LAMBDA_B__LAMBDA_L_L)

        .value("PHI_D", Observables::PHI_D)
        .value("DELTA_M_BD", Observables::DELTA_M_BD)

        .value("PHI_S", Observables::PHI_S)
        .value("DELTA_M_BS", Observables::DELTA_M_BS)
        .value("A_FS", Observables::A_FS)

        .value("DELTA_M_K", Observables::DELTA_M_K)
        .value("ABS_EPSILON_K", Observables::ABS_EPSILON_K)

        .value("X_D", Observables::X_D)

        .value("BR_KL__MU_MU", Observables::BR_KL__MU_MU)
        .value("BR_KS__MU_MU", Observables::BR_KS__MU_MU)

        .value("BR_K__MU_NU__BR_PI__MU_NU", Observables::BR_K__MU_NU__BR_PI__MU_NU)
        .value("R_MU23", Observables::R_MU23)

        .value("BR_K__PI_NU_NU", Observables::BR_K__PI_NU_NU)
        .value("BR_KL__PI0_NU_NU", Observables::BR_KL__PI0_NU_NU)

        .value("BR_D__MU_NU", Observables::BR_D__MU_NU)
        .value("BR_DS__MU_NU", Observables::BR_DS__MU_NU)
        .value("BR_DS__TAU_NU", Observables::BR_DS__TAU_NU)


        .export_values();

    py::enum_<Decays>(m, "Decays")
        .value("B__D_l_nu", Decays::B__D_l_nu)
        .value("B__Dstar_l_nu", Decays::B__Dstar_l_nu)
        .value("B__Kstar", Decays::B__Kstar_gamma)
        .value("B__l_l", Decays::B__l_l)
        .value("B__l_nu", Decays::B__l_nu)
        .value("B__Xs_gamma", Decays::B__Xs_gamma)
        .value("B__Xs_l_l", Decays::B__Xs_l_l)
        .value("B__K_l_l", Decays::B__K_l_l)
        .value("B__Kstar_l_l", Decays::B__Kstar_l_l)
        .value("Bs__phi_l_l", Decays::Bs__phi_l_l)
        .value("Lambda_b__Lambda_l_l", Decays::Lambda_b__Lambda_l_l)
        .value("M0_Mix", Decays::M0_Mix)
        .value("K__l_l", Decays::K__l_l)
        .value("K__pi_nu_nu", Decays::K__pi_nu_nu)
        .value("K__l_nu", Decays::K__l_nu)
        .value("D__l_nu", Decays::D__l_nu)
        .value("Ds__l_nu", Decays::Ds__l_nu)
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
        .value("C_V1_bc", WCoef::C_V1_bc)
        .value("C_V2_bc", WCoef::C_V2_bc)
        .value("C_S1_bc", WCoef::C_S1_bc)
        .value("C_S2_bc", WCoef::C_S2_bc)
        .value("C_T_bc", WCoef::C_T_bc)
        .value("C_V1_bu", WCoef::C_V1_bu)
        .value("C_V2_bu", WCoef::C_V2_bu)
        .value("C_S1_bu", WCoef::C_S1_bu)
        .value("C_S2_bu", WCoef::C_S2_bu)
        .value("C_T_bu", WCoef::C_T_bu)
        .value("C_V1_cs", WCoef::C_V1_cs)
        .value("C_V2_cs", WCoef::C_V2_cs)
        .value("C_S1_cs", WCoef::C_S1_cs)
        .value("C_S2_cs", WCoef::C_S2_cs)
        .value("C_T_cs", WCoef::C_T_cs)
        .value("C_V1_cd", WCoef::C_V1_cd)
        .value("C_V2_cd", WCoef::C_V2_cd)
        .value("C_S1_cd", WCoef::C_S1_cd)
        .value("C_S2_cd", WCoef::C_S2_cd)
        .value("C_T_cd", WCoef::C_T_cd)
        .value("C_V1_su", WCoef::C_V1_su)
        .value("C_V2_su", WCoef::C_V2_su)
        .value("C_S1_su", WCoef::C_S1_su)
        .value("C_S2_su", WCoef::C_S2_su)
        .value("C_T_su", WCoef::C_T_su)
        .value("C_V1_du", WCoef::C_V1_du)
        .value("C_V2_du", WCoef::C_V2_du)
        .value("C_S1_du", WCoef::C_S1_du)
        .value("C_S2_du", WCoef::C_S2_du)
        .value("C_T_du", WCoef::C_T_du)
        .value("C_BD_1", WCoef::C_BD_1)
        .value("CT_BD_1", WCoef::CT_BD_1)
        .value("C_BD_2", WCoef::C_BD_2)
        .value("CT_BD_2", WCoef::CT_BD_2)
        .value("C_BD_3", WCoef::C_BD_3)
        .value("CT_BD_3", WCoef::CT_BD_3)
        .value("C_BD_4", WCoef::C_BD_4)
        .value("C_BD_5", WCoef::C_BD_5)
        .value("C_BS_1", WCoef::C_BS_1)
        .value("CT_BS_1", WCoef::CT_BS_1)
        .value("C_BS_2", WCoef::C_BS_2)
        .value("CT_BS_2", WCoef::CT_BS_2)
        .value("C_BS_3", WCoef::C_BS_3)
        .value("CT_BS_3", WCoef::CT_BS_3)
        .value("C_BS_4", WCoef::C_BS_4)
        .value("C_BS_5", WCoef::C_BS_5)
        .value("C_SD_1", WCoef::C_SD_1)
        .value("CT_SD_1", WCoef::CT_SD_1)
        .value("C_SD_2", WCoef::C_SD_2)
        .value("CT_SD_2", WCoef::CT_SD_2)
        .value("C_SD_3", WCoef::C_SD_3)
        .value("CT_SD_3", WCoef::CT_SD_3)
        .value("C_SD_4", WCoef::C_SD_4)
        .value("C_SD_5", WCoef::C_SD_5)
        .value("C_CU_1", WCoef::C_CU_1)
        .value("CT_CU_1", WCoef::CT_CU_1)
        .value("C_CU_2", WCoef::C_CU_2)
        .value("CT_CU_2", WCoef::CT_CU_2)
        .value("C_CU_3", WCoef::C_CU_3)
        .value("CT_CU_3", WCoef::CT_CU_3)
        .value("C_CU_4", WCoef::C_CU_4)
        .value("C_CU_5", WCoef::C_CU_5)
        .value("CK9", WCoef::CK9)
        .value("CPK9", WCoef::CPK9)
        .value("CK10", WCoef::CK10)
        .value("CPK10", WCoef::CPK10)
        .value("CKQ1", WCoef::CKQ1)
        .value("CKQ2", WCoef::CKQ2)
        .value("CPKQ1", WCoef::CPKQ1)
        .value("CPKQ2", WCoef::CPKQ2)
        .value("CK_L", WCoef::CK_L)
        .export_values();

    py::enum_<WGroup>(m, "WGroup")
        .value("B", WGroup::B)
        .value("BPrime", WGroup::BPrime)
        .value("BScalar", WGroup::BScalar)
        .value("CC_bc", WGroup::CC_bc)
        .value("CC_bu", WGroup::CC_bu)
        .value("CC_cs", WGroup::CC_cs)
        .value("CC_cd", WGroup::CC_cd)
        .value("CC_su", WGroup::CC_su)
        .value("CC_du", WGroup::CC_du)
        .value("MESON_MIXING", WGroup::MESON_MIXING)
        .value("K", WGroup::K)
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

    py::enum_<DataType>(m, "DataType")
        .value("VALUE", DataType::VALUE)
        .value("STD_STAT", DataType::STD_STAT)
        .value("STD_SYST", DataType::STD_SYST)
        .value("STD_COMBINED", DataType::STD_COMBINED)
        .export_values();

    py::enum_<UncertaintyType>(m, "UncertaintyType")
        .value("STAT", UncertaintyType::STAT)
        .value("SYST", UncertaintyType::SYST)
        .value("COMBINED", UncertaintyType::COMBINED)
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
    
    // BIND_ENUM_MAPPER(ObservableMapper, Observables)
    BIND_ENUM_MAPPER(OrderMapper, QCDOrder)
    // BIND_ENUM_MAPPER(GroupMapper, WGroup)
    // BIND_ENUM_MAPPER(WCoefMapper, WCoef)
    BIND_ENUM_MAPPER(ParameterTypeMapper, ParameterType)
    BIND_ENUM_MAPPER(ModelMapper, Model)
    BIND_ENUM_MAPPER(WilsonBasisMapper, WilsonBasis)
    BIND_ENUM_MAPPER(ContributionTypeMapper, ContributionType)
    BIND_ENUM_MAPPER(MassTypeMapper, MassType)
    // BIND_ENUM_MAPPER(ScaleTypeMapper, ScaleType)
    // BIND_ENUM_MAPPER(DecayMapper, Decays)

    py::class_<ObservableId>(m, "ObservableId")
        .def(py::init<>())
        .def("__str__", &ObservableId::str)
        .def("str", &ObservableId::str);

    py::class_<ObservableMapper, std::shared_ptr<ObservableMapper>>(m, "ObservableMapper")
        .def_static("str",
            py::overload_cast<Observables>(&ObservableMapper::str),
            py::arg("obs"))
        .def_static("enum_elt",   &ObservableMapper::enum_elt_legacy, py::arg("name"))
        .def_static("get_str",    &ObservableMapper::get_str)
        .def_static("get_enum",   &ObservableMapper::get_enum)
        .def_static("get_str_all",&ObservableMapper::get_str_all)

        // runtime id/ext
        .def_static("id_of",        &ObservableMapper::id_of,         py::arg("name"))
        .def_static("canonical",
            py::overload_cast<const ObservableId&>(&ObservableMapper::str),
            py::arg("id"))
        .def_static("from_flha",    &ObservableMapper::from_flha,     py::arg("lha"))

        // ⬇️ DISAMBIGUATION: bind the "id -> LhaID" overload
        .def_static("flha",
            py::overload_cast<const ObservableId&>(&ObservableMapper::flha),
            py::arg("id"));

    // py::class_<WCoefMapper, std::shared_ptr<WCoefMapper>>(m, "WCoefMapper")
    //     .def_static("str", &WCoefMapper::str, py::arg("coef"))
    //     .def_static("enum_elt", &WCoefMapper::enum_elt, py::arg("name"))
    //     .def_static("get_str", &WCoefMapper::get_str)
    //     .def_static("get_enum", &WCoefMapper::get_enum)
    //     .def_static("get_group", &WCoefMapper::get_group, py::arg("group"))
    //     .def_static("flha_base", &WCoefMapper::flha_base, py::arg("coef"))
    //     .def_static("flha_full", &WCoefMapper::flha_full, py::arg("coef"), py::arg("order"), py::arg("type"))
    //     .def_static("from_flha", &WCoefMapper::from_flha, py::arg("content"), py::arg("structure"))
    //     .def_static("n_wilsons", &WCoefMapper::n_wilsons)
    //     .def_static("mapping", &WCoefMapper::mapping, py::return_value_policy::reference)
    //     .def_static("inverse_mapping", &WCoefMapper::inverse_mapping, py::return_value_policy::reference)
    //     .def_static("flha_mapping", &WCoefMapper::flha_mapping, py::return_value_policy::reference)
    //     .def_static("inverse_flha_mapping", &WCoefMapper::inverse_flha_mapping, py::return_value_policy::reference)
    //     .def_static("B_group", &WCoefMapper::B_group, py::return_value_policy::reference)
    //     .def_static("B_prime_group", &WCoefMapper::B_prime_group, py::return_value_policy::reference)
    //     .def_static("B_scalar_group", &WCoefMapper::B_scalar_group, py::return_value_policy::reference)
    //     // .def_static("B_lnu_group", &WCoefMapper::B_lnu_group, py::return_value_policy::reference)
    //     .def_static("b_clnu_group", &WCoefMapper::b_clnu_group, py::return_value_policy::reference);

    py::class_<WCoefMapper, std::shared_ptr<WCoefMapper>>(m, "WCoefMapper")
        // legacy: str(enum) + enum_elt(name)->Enum
        .def_static("str",
            py::overload_cast<WCoef>(&WCoefMapper::str),
            py::arg("coef"))
        .def_static("enum_elt",   &WCoefMapper::enum_elt_legacy, py::arg("name"))
        .def_static("get_str",    &WCoefMapper::get_str)
        .def_static("get_enum",   &WCoefMapper::get_enum)
        .def_static("get_str_all",&WCoefMapper::get_str_all)

        // runtime id API
        .def_static("id_of",        &WCoefMapper::id_of,        py::arg("name"))
        .def_static("canonical",
            py::overload_cast<const WCoefId&>(&WCoefMapper::str),
            py::arg("id"))
        .def_static("from_external",&WCoefMapper::from_external,py::arg("flha_pair"))
        .def_static("external_of",  &WCoefMapper::external_of,  py::arg("id"))
        .def_static("register_custom",
            [](const std::string& canon, const std::vector<std::string>& aliases, std::pair<int,int> flha){
                return WCoefMapper::register_custom(canon, aliases, flha);
            },
            py::arg("canonical"), py::arg("aliases") = std::vector<std::string>{}, py::arg("flha"))

        // groupes legacy
        .def_static("get_group",        &WCoefMapper::get_group,        py::arg("group"))
        .def_static("B_group",          &WCoefMapper::B_group,          py::return_value_policy::reference)
        .def_static("B_prime_group",    &WCoefMapper::B_prime_group,    py::return_value_policy::reference)
        .def_static("B_scalar_group",   &WCoefMapper::B_scalar_group,   py::return_value_policy::reference)
        .def_static("b_clnu_group",     &WCoefMapper::b_clnu_group,     py::return_value_policy::reference)

        // FLHA legacy — ⬇️ DEUX overloads, on les disambiguë
        .def_static("flha_base",
            py::overload_cast<WCoef>(&WCoefMapper::flha_base),
            py::arg("coef"))
        .def_static("flha_base",
            py::overload_cast<const WCoefId&>(&WCoefMapper::flha_base),
            py::arg("id"));



    // py::class_<GroupMapper, std::shared_ptr<GroupMapper>>(m, "GroupMapper")
    //     .def_static(
    //         "str",
    //         static_cast<std::string(*)(WGroup)>(&GroupMapper::str),
    //         py::arg("group")
    //     )
    //     .def_static(
    //         "str",
    //         static_cast<std::string(*)(WGroup, ScaleType, WilsonBasis)>(&GroupMapper::str),
    //         py::arg("group"), py::arg("scale"), py::arg("basis") = WilsonBasis::B_STANDARD
    //     )
    //     .def_static("enum_elt", &GroupMapper::enum_elt, py::arg("name"))
    //     .def_static("get_str", &GroupMapper::get_str)
    //     .def_static("get_enum", &GroupMapper::get_enum);

    py::class_<WGroupId>(m, "WGroupId")
        .def(py::init<>())
        // pour que str(id) marche en Python
        .def("__str__", &WGroupId::str)
        .def("str", &WGroupId::str);

    py::class_<GroupMapper, std::shared_ptr<GroupMapper>>(m, "GroupMapper")
        // str(enum) et str(enum, scale, basis)
        .def_static("str", py::overload_cast<WGroup>(&GroupMapper::str), py::arg("group"))
        .def_static("str", py::overload_cast<WGroup, ScaleType, WilsonBasis>(&GroupMapper::str),
                    py::arg("group"), py::arg("scale"), py::arg("basis") = WilsonBasis::B_STANDARD)

        // legacy enum_elt(name)->WGroup + listes
        .def_static("enum_elt",     &GroupMapper::enum_elt_legacy,  py::arg("name"))
        .def_static("get_str",      &GroupMapper::get_str)
        .def_static("get_enum",     &GroupMapper::get_enum)
        .def_static("get_str_all",  &GroupMapper::get_str_all)

        // (optionnel) runtime id API
        .def_static("id_of",     &GroupMapper::id_of,     py::arg("name"))
        .def_static("canonical", py::overload_cast<const WGroupId&>(&GroupMapper::str), py::arg("id"));


    // py::class_<ScaleTypeMapper, std::shared_ptr<ScaleTypeMapper>>(m, "ScaleTypeMapper")
    //     .def_static("str", &ScaleTypeMapper::str, py::arg("type"))
    //     .def_static("enum_elt", &ScaleTypeMapper::enum_elt, py::arg("name"))
    //     .def_static("get_str", &ScaleTypeMapper::get_str)
    //     .def_static("get_enum", &ScaleTypeMapper::get_enum)
    //     .def_static("block", &ScaleTypeMapper::block, py::arg("type"));

    py::class_<ScaleTypeMapper, std::shared_ptr<ScaleTypeMapper>>(m, "ScaleTypeMapper")
        .def_static("str",        py::overload_cast<ScaleType>(&ScaleTypeMapper::str), py::arg("type"))
        .def_static("enum_elt",   &ScaleTypeMapper::enum_elt_legacy, py::arg("name"))
        .def_static("get_str",    &ScaleTypeMapper::get_str)
        .def_static("get_enum",   &ScaleTypeMapper::get_enum)
        .def_static("get_str_all",&ScaleTypeMapper::get_str_all)
        .def_static("block",      &ScaleTypeMapper::block, py::arg("type"))
        // (optionnel) runtime id:
        .def_static("id_of",      &ScaleTypeMapper::id_of, py::arg("name"))
        .def_static("canonical",  py::overload_cast<const ScaleTypeId&>(&ScaleTypeMapper::str), py::arg("id"));

    // py::class_<DecayMapper, std::shared_ptr<DecayMapper>>(m, "DecayMapper")
    //     .def_static("str", &DecayMapper::str, py::arg("type"))
    //     .def_static("enum_elt", &DecayMapper::enum_elt, py::arg("name"))
    //     .def_static("get_str", &DecayMapper::get_str)
    //     .def_static("get_enum", &DecayMapper::get_enum)
    //     .def_static("get_observables", &DecayMapper::get_observables, py::arg("decay"))
    //     .def_static("get_decay", &DecayMapper::get_decay, py::arg("observable"));

    py::class_<DecayMapper, std::shared_ptr<DecayMapper>>(m, "DecayMapper")
        .def_static("str",        py::overload_cast<Decays>(&DecayMapper::str), py::arg("type"))
        .def_static("enum_elt",   &DecayMapper::enum_elt_legacy, py::arg("name"))
        .def_static("get_str",    &DecayMapper::get_str)
        .def_static("get_enum",   &DecayMapper::get_enum)
        .def_static("get_str_all",&DecayMapper::get_str_all)

        // métier (comme avant)
        .def_static("get_observables", &DecayMapper::get_observables, py::arg("decay"))
        .def_static("get_decay",       &DecayMapper::get_decay,       py::arg("observable"))

        // (optionnel) runtime id/ext
        .def_static("id_of",        &DecayMapper::id_of,         py::arg("name"))
        .def_static("canonical",    py::overload_cast<const DecayId&>(&DecayMapper::str), py::arg("id"))
        .def_static("from_external",&DecayMapper::from_external, py::arg("lha"))
        .def_static("external_of",  &DecayMapper::external_of,   py::arg("id"))
        .def_static("set_external", &DecayMapper::set_external,  py::arg("id"), py::arg("lha"));
 
    py::class_<LhaID>(m, "LhaID")
    .def(py::init<const std::string&>(), py::arg("parts"))
    .def(py::init<long>(), py::arg("id"))
    .def(py::init<const std::vector<long>&>(), py::arg("sub_ids"))
    .def(py::init<std::initializer_list<long>>(), py::arg("sub_ids"))
    .def(py::init([](py::args args) {
        std::vector<long> values;
        for (auto& item : args) {
            values.push_back(item.cast<long>());
        }
        return LhaID(values);
    }), "Construct from multiple long arguments")

    .def("to_string", &LhaID::to_string)
    .def("get_parts", &LhaID::get_parts)

    .def("__int__", [](const LhaID& self) { return static_cast<long>(self); })

    .def("__repr__", [](const LhaID& self) {
        return "<LhaID: " + self.to_string() + ">";
    })

    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self < py::self)

    .def("__hash__", [](const LhaID& self) {
        return std::hash<LhaID>{}(self);
    });

    py::class_<BlockName>(m, "BlockName")

    .def(py::init<>())
    .def(py::init<const std::string&>(), py::arg("name"))
    .def(py::init<const char*>(), py::arg("name"))
    .def(py::init<std::initializer_list<std::string>>(), py::arg("names"))
    .def(py::init<const std::unordered_set<std::string>&>(), py::arg("names"))

    .def("get_alias", &BlockName::get_alias)
    .def("has_alias", &BlockName::hasAlias, py::arg("alias"))
    .def("add_alias", &BlockName::addAlias, py::arg("alias"), py::return_value_policy::reference)
    .def("to_string", &BlockName::to_string)
    .def("to_upper", &BlockName::to_upper)

    .def("__str__", &BlockName::to_string)
    .def("__repr__", [](const BlockName& b) {
        std::ostringstream oss;
        oss << "<BlockName: ";
        bool first = true;
        for (const auto& name : b.get_alias()) {
            if (!first) oss << "/";
            oss << name;
            first = false;
        }
        oss << ">";
        return oss.str();
    })

    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self < py::self)

    .def("__eq__", [](const BlockName& self, const std::string& s) { return self == s; })
    .def("__ne__", [](const BlockName& self, const std::string& s) { return self != s; })
    .def("__eq__", [](const std::string& s, const BlockName& self) { return self == s; })
    .def("__ne__", [](const std::string& s, const BlockName& self) { return self != s; })

    .def("__hash__", [](const BlockName& b) {
        return std::hash<BlockName>{}(b);
    });

    py::class_<ParamId>(m, "ParamId")
    .def(py::init<>())
    .def(py::init<const BlockName&, const LhaID&>(), py::arg("block"), py::arg("code"))
    .def(py::init<ParameterType, const BlockName&, const LhaID&>(), py::arg("type"), py::arg("block"), py::arg("code"))

    .def_readwrite("type", &ParamId::type)
    .def_readwrite("block", &ParamId::block)
    .def_readwrite("code", &ParamId::code)

    .def("set_parameter_type", &ParamId::set_parameter_type)

    .def(py::self == py::self)
    .def(py::self < py::self)

    .def("__repr__", [](const ParamId& pid) {
        std::ostringstream oss;
        oss << "<ParamId: " << pid.block << ":" << pid.code;
        if (pid.type.has_value()) {
            oss << ", type=" << static_cast<int>(pid.type.value()); // adapt if enum is bound
        } else {
            oss << ", type=None";
        }
        oss << ">";
        return oss.str();
    })

    .def("__hash__", [](const ParamId& pid) {
        return std::hash<ParamId>{}(pid);
    });

    // py::class_<DependenciesHelper, std::shared_ptr<DependenciesHelper>>(m, "DependenciesHelper")
    //     .def_static("get_allowed_parameters", &DependenciesHelper::get_allowed_parameters, py::arg("obs"))
    //     .def_static("is_param_allowed", &DependenciesHelper::is_param_allowed, py::arg("obs"), py::arg("param"));

    py::class_<DependenciesHelper, std::shared_ptr<DH>>(m, "DependenciesHelper")
    // version avec Observables
    .def_static(
        "get_allowed_parameters",
        py::overload_cast<Observables>(&DH::get_allowed_parameters),
        py::arg("obs")
    )
    // autre nom python pour la version ObservableId
    .def_static(
        "get_allowed_parameters_from_id",
        py::overload_cast<ObservableId>(&DH::get_allowed_parameters),
        py::arg("obs_id")
    )
    // version avec Observables
    .def_static(
        "is_param_allowed",
        py::overload_cast<Observables, ParamId>(&DH::is_param_allowed),
        py::arg("obs"),
        py::arg("param")
    )
    // autre nom pour ObservableId
    .def_static(
        "is_param_allowed_from_id",
        py::overload_cast<ObservableId, ParamId>(&DH::is_param_allowed),
        py::arg("obs_id"),
        py::arg("param")
    );

    py::class_<LhaParamsHelper, std::shared_ptr<LhaParamsHelper>>(m, "LhaParamsHelper")
        .def_static("get_minimal_content", &LhaParamsHelper::get_minimal_content, py::arg("block_name"));


    py::class_<WilsonBuildConfig>(m, "WilsonBuildConfig")
        .def(py::init<>())
        .def_readwrite("groups", &WilsonBuildConfig::groups)
        .def_readwrite("matching_scale", &WilsonBuildConfig::matching_scale)
        .def_readwrite("hadronic_scale", &WilsonBuildConfig::hadronic_scale)
        .def_readwrite("order", &WilsonBuildConfig::order);
    
    py::class_<WilsonRequest>(m, "WilsonRequest")
        .def(py::init<WGroup, WCoef, QCDOrder, ContributionType, ScaleType, bool>())
        .def_readwrite("group", &WilsonRequest::group)
        .def_readwrite("coefficient", &WilsonRequest::coefficient)
        .def_readwrite("order", &WilsonRequest::order)
        .def_readwrite("contribution", &WilsonRequest::contribution)
        .def_readwrite("scale_type", &WilsonRequest::scale_type)
        .def_readwrite("sum_qcd_orders", &WilsonRequest::sum_qcd_orders);
    
    py::class_<AbstractConfig>(m, "AbstractConfig"); 

    py::class_<AlphasConfig, AbstractConfig>(m, "AlphasConfig")
        .def(py::init<double, MassType, MassType>(), py::arg("scale"), py::arg("m_b_type"), py::arg("m_t_type"))
        .def_readwrite("scale", &AlphasConfig::scale)
        .def_readwrite("m_b_type", &AlphasConfig::m_b_type)
        .def_readwrite("m_t_type", &AlphasConfig::m_t_type);

    py::class_<MassConfig, AlphasConfig>(m, "MassConfig")
        .def(py::init<int, double, MassType, MassType>(), py::arg("pdg_id"), py::arg("scale"), py::arg("m_b_type"), py::arg("m_t_type"))
        .def_readwrite("pdg_id", &MassConfig::pdg_id);

}
