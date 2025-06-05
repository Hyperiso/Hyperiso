

#include "Map.h"

const std::map<Observables, std::string>& observable_mapping() {
    static const std::map<Observables, std::string> m = {
        {Observables::BR_BS_MUMU, "BR_Bs__mu_mu"},
        {Observables::BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
        {Observables::BR_BD_MUMU, "BR_Bd__mu_mu"},
        {Observables::BR_BU_TAU_NU, "BR_Bu__tau_nu"},
        {Observables::R_TAU_NU, "R_tau_nu"},
        {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B__K*_gamma"},
        {Observables::BR_B_XS_GAMMA, "BR_B__Xs_gamma"},
        {Observables::BR_B__D_TAU_NU, "BR_B__D_tau_nu"},
        {Observables::A_FB_B__D_TAU_NU, "A_FB_B__D_tau_nu"},
        {Observables::P_TAU_B__D_TAU_NU, "P_tau_B__D_tau_nu"},
        {Observables::R_D, "R_D"},
        {Observables::BR_B__DSTAR_TAU_NU, "BR_B__D*_tau_nu"},
        {Observables::A_FB_B__DSTAR_TAU_NU, "A_FB_B__D*_tau_nu"},
        {Observables::P_TAU_B__DSTAR_TAU_NU, "P_tau_B__D*_tau_nu"},
        {Observables::P_D_B__DSTAR_TAU_NU, "P_D*_B__D*_tau_nu"},
        {Observables::R_DSTAR, "R_D*"},
    };
    return m;
}

const std::map<QCDOrder, std::string>& order_mapping() {
    static const std::map<QCDOrder, std::string> m = {
        {QCDOrder::NONE, "None"},
        {QCDOrder::LO, "LO"},
        {QCDOrder::NLO, "NLO"},
        {QCDOrder::NNLO, "NNLO"},
    };
    return m;
}

const std::map<WGroup, std::string>& group_mapping() {
    static const std::map<WGroup, std::string> m = {
        {WGroup::B, "BCoefficients"},
        {WGroup::BPrime, "BPrimeCoefficients"},
        {WGroup::BScalar, "BScalarCoefficients"},
        {WGroup::BCC, "BChargedCurrentCoefficients"},
    };
    return m;
}

const std::map<WCoef, std::string>& wcoef_mapping() {
    static const std::map<WCoef, std::string> m = {
        {WCoef::C1, "C1"}, {WCoef::C2, "C2"}, {WCoef::C3, "C3"}, {WCoef::C4, "C4"}, {WCoef::C5, "C5"}, {WCoef::C6, "C6"},
        {WCoef::C7, "C7"}, {WCoef::C8, "C8"}, {WCoef::C9, "C9"}, {WCoef::C10, "C10"},
        {WCoef::CP1, "CP1"}, {WCoef::CP2, "CP2"}, {WCoef::CP3, "CP3"}, {WCoef::CP4, "CP4"},
        {WCoef::CP5, "CP5"}, {WCoef::CP6, "CP6"}, {WCoef::CP7, "CP7"}, {WCoef::CP8, "CP8"},
        {WCoef::CP9, "CP9"}, {WCoef::CP10, "CP10"}, {WCoef::CQ1, "CQ1"}, {WCoef::CQ2, "CQ2"},
        {WCoef::CPQ1, "CPQ1"}, {WCoef::CPQ2, "CPQ2"},
        {WCoef::C_V1, "C_V1"}, {WCoef::C_V2, "C_V2"}, {WCoef::C_S1, "C_S1"}, {WCoef::C_S2, "C_S2"}, {WCoef::C_T, "C_T"},
    };
    return m;
}

const std::map<WCoef, std::pair<int, int>>& wcoef_flha_mapping() {
    static const std::map<WCoef, std::pair<int, int>> m = {
        {WCoef::C1, {3040405, 6161}}, {WCoef::C2, {3040405, 4141}}, {WCoef::C3, {3050707, 4133}},
        {WCoef::C4, {3050707, 6153}}, {WCoef::C5, {3050707, 4536}}, {WCoef::C6, {3050707, 6556}},
        {WCoef::C7, {305, 4422}}, {WCoef::C8, {305, 6421}}, {WCoef::C9, {3051313, 4133}}, {WCoef::C10, {3051313, 4137}},
        {WCoef::CP1, {3040405, 6262}}, {WCoef::CP2, {3040405, 4242}}, {WCoef::CP3, {3050707, 4233}},
        {WCoef::CP4, {3050707, 6253}}, {WCoef::CP5, {3050707, 4636}}, {WCoef::CP6, {3050707, 6656}},
        {WCoef::CP7, {305, 4322}}, {WCoef::CP8, {305, 4321}}, {WCoef::CP9, {3051313, 4233}},
        {WCoef::CP10, {3051313, 4234}}, {WCoef::CQ1, {3051313, 3230}}, {WCoef::CQ2, {3051313, 3233}},
        {WCoef::CPQ1, {3051313, 3130}}, {WCoef::CPQ2, {3051313, 3133}},
        {WCoef::C_V1, {4051516, 4141}}, {WCoef::C_V2, {4051516, 4241}},
        {WCoef::C_S1, {4051516, 3231}}, {WCoef::C_S2, {4051516, 3131}}, {WCoef::C_T, {4051516, 4343}},
    };
    return m;
}


const std::map<ParameterType, std::string>& parametertype_mapping() {
    static const std::map<ParameterType, std::string> m = {
        {ParameterType::SM, "SM"},
        {ParameterType::BSM, "BSM"},
        {ParameterType::FLAVOR, "FLAVOR"},
        {ParameterType::WILSON, "WILSON"},
        {ParameterType::DECAY, "DECAY"},
        {ParameterType::PASSTHROUGH, "PASSTHROUGH"},
        {ParameterType::OBSERVABLE, "OBSERVABLE"},
    };
    return m;
}

const std::map<Model, std::string>& model_mapping() {
    static const std::map<Model, std::string> m = {
        {Model::SM, "SM"},
        {Model::SUSY, "SUSY"},
        {Model::THDM, "THDM"},
        {Model::CUSTOM, "CUSTOM"},
    };
    return m;
}

const std::map<WilsonBasis, std::string>& wilsonbasis_mapping() {
    static const std::map<WilsonBasis, std::string> m = {
        {WilsonBasis::B_STANDARD, "STANDARD"},
        {WilsonBasis::B_TRADITIONAL, "TRADITIONAL"},
    };
    return m;
}

const std::map<ContributionType, std::string>& contributiontype_mapping() {
    static const std::map<ContributionType, std::string> m = {
        {ContributionType::SM, "SM"},
        {ContributionType::BSM, "BSM"},
        {ContributionType::TOTAL, "TOTAL"},
    };
    return m;
}

const std::map<Decays, std::string>& decays_mapping() {
    static const std::map<Decays, std::string> m = {
        {Decays::B__D_l_nu, "B__D_l_nu"},
        {Decays::B__Dstar_l_nu, "B__Dstar_l_nu"},
        {Decays::B__Kstar, "B__Kstar"},
        {Decays::B__l_l, "B__l_l"},
        {Decays::B__l_nu, "B__l_nu"},
        {Decays::B__Xs, "B__Xs"},
    };
    return m;
}

const std::map<Decays, std::vector<Observables>>& decay_observable_mapping() {
    static const std::map<Decays, std::vector<Observables>> m = {
        {Decays::B__D_l_nu,     {Observables::BR_B__D_TAU_NU, Observables::A_FB_B__D_TAU_NU, Observables::P_TAU_B__D_TAU_NU, Observables::R_D}},
        {Decays::B__Dstar_l_nu, {Observables::BR_B__DSTAR_TAU_NU, Observables::A_FB_B__DSTAR_TAU_NU, Observables::P_TAU_B__DSTAR_TAU_NU, Observables::P_D_B__DSTAR_TAU_NU, Observables::R_DSTAR}},
        {Decays::B__Kstar,      {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA}},
        {Decays::B__l_l,        {Observables::BR_BS_MUMU, Observables::BR_BD_MUMU, Observables::BR_BS_MUMU_UNTAG}},
        {Decays::B__l_nu,       {Observables::BR_BU_TAU_NU, Observables::R_TAU_NU}},
        {Decays::B__Xs,         {Observables::BR_B_XS_GAMMA}},
    };
    return m;
}

