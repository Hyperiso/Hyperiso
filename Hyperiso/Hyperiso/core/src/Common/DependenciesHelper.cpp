#include "DependenciesHelper.h"


const std::map<DecayId, std::unordered_set<ParamId>> DependenciesHelper::dep_lists = {
    {DecayMapper::to_id(Decays::B__D_l_nu), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "MASS", 11},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "VCKM", {1, 2}},
        ParamId{ParameterType::FLAVOR, "FMASS", 511},
        ParamId{ParameterType::FLAVOR, "FMASS", 411},
        ParamId{ParameterType::FLAVOR, "FLIFE", 511},
        ParamId{ParameterType::DECAY, "B_Dlnu", 1},
        ParamId{ParameterType::DECAY, "B_Dlnu", 2},
        ParamId{ParameterType::DECAY, "B_Dlnu", 3},
        // TODO : Wilson all CC_bc group
        // TODO : What about QCD ?
    }},
    {DecayMapper::to_id(Decays::B__Dstar_l_nu), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "MASS", 11},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "VCKM", {1, 2}},
        ParamId{ParameterType::FLAVOR, "FMASS", 511},
        ParamId{ParameterType::FLAVOR, "FMASS", 413},
        ParamId{ParameterType::FLAVOR, "FLIFE", 511},
        ParamId{ParameterType::DECAY, "B_Dslnu", 1},
        ParamId{ParameterType::DECAY, "B_Dslnu", 2},
        ParamId{ParameterType::DECAY, "B_Dslnu", 3},
        ParamId{ParameterType::DECAY, "B_Dslnu", 4},
        // TODO : Wilson all CC_bc group
    }},
    {DecayMapper::to_id(Decays::B__K_l_l), {
        
    }},
    {DecayMapper::to_id(Decays::B__Kstar_l_l), {
        ParamId{ParameterType::SM, "EW", {1, 2}},
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "MASS", 11},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "QCD", {5, 2}},
        ParamId{ParameterType::SM, "VCKM", {0, 1}},
        ParamId{ParameterType::SM, "VCKM", {0, 2}},
        ParamId{ParameterType::SM, "VCKM", {2, 1}},
        ParamId{ParameterType::SM, "VCKM", {2, 2}},
        ParamId{ParameterType::FLAVOR, "FMASS", 511},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FLIFE", 511},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 313},
        ParamId{ParameterType::FLAVOR, "FMASS", 323},
        ParamId{ParameterType::FLAVOR, "FCONST", {313, 1}},
        ParamId{ParameterType::FLAVOR, "FCONST", {313, 2}},
        ParamId{ParameterType::FLAVOR, "FCONST", {323, 1}},
        ParamId{ParameterType::FLAVOR, "FCONST", {323, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {7, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {7, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {8, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {8, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", 9},
        ParamId{ParameterType::DECAY, "B_Ks", 10},
        ParamId{ParameterType::DECAY, "B_Ks", 11},
        ParamId{ParameterType::DECAY, "B_Ks", {12, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {12, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", 14},
        ParamId{ParameterType::DECAY, "B_Ks", {15, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {15, 2}},
        // ParamId{ParameterType::DECAY, "B_Ks", {1, 0, 1}},
        // ParamId{ParameterType::DECAY, "B_Ks", {1, 0, 2}},
        // ParamId{ParameterType::DECAY, "B_Ks", {1, 0, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 1, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 1, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 1, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 2, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 2, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 2, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 3, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 3, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 3, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 4, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 4, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 4, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 5, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 5, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 5, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 6, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 6, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 6, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 7, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 7, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {1, 7, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 1, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 1, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 1, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 2, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 2, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 2, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 3, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 3, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 4, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 4, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 4, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 5, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 5, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 5, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 6, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 6, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 7, 0}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 7, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {3, 7, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 1, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 1, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 1, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 1, 4}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 1, 5}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 1, 6}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 2, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 2, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 2, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 2, 4}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 2, 5}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 2, 6}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 4}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 5}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 6}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 7}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 3, 8}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 4, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 4, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 4, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 4, 4}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 4, 5}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 4, 6}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 5, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 5, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 5, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 5, 4}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 5, 5}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 5, 6}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 1}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 2}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 3}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 4}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 5}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 6}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 7}},
        ParamId{ParameterType::DECAY, "B_Ks", {18, 6, 8}},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C1, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C2, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C3, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C4, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C5, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C6, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C7, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C8, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C1, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C2, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C3, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C4, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C5, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C6, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C7, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C8, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C9, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C10, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C1, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C2, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C3, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C4, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C5, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C6, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C7, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C8, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C9, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::B, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::C10, QCDOrder::NNLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CP7, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CP9, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CP10, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CPQ1, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CPQ2, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CQ1, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CQ2, QCDOrder::LO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CQ1, QCDOrder::NLO, ContributionType::SM)},
        ParamId{ParameterType::WILSON, GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC), WCoefMapper::flha_full(WCoef::CQ2, QCDOrder::NLO, ContributionType::SM)},
    }},
    {DecayMapper::to_id(Decays::B__Kstar_gamma), {
        
    }},
    {DecayMapper::to_id(Decays::B__l_l), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "EW", {1, 2}},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "VCKM", {2, 0}},
        ParamId{ParameterType::SM, "VCKM", {2, 1}},
        ParamId{ParameterType::SM, "VCKM", {2, 2}},
        ParamId{ParameterType::SM, "QCD", {5, 2}},
        ParamId{ParameterType::FLAVOR, "FMASS", 531},
        ParamId{ParameterType::FLAVOR, "FMASS", 511},
        ParamId{ParameterType::FLAVOR, "FLIFE", 531},
        ParamId{ParameterType::FLAVOR, "FLIFE", 511},
        ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}},
        ParamId{ParameterType::FLAVOR, "FCONST", {531, 1}},
        ParamId{ParameterType::DECAY, "B_ll", 1},
        ParamId{ParameterType::DECAY, "B_ll", 2},
        // ParamId{ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD", {3051313, 4137, 0, 2}},
        // ParamId{ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD", {3051313, 4137, 1, 2}},
        // ParamId{ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD", {3051313, 4137, 2, 2}},
        // ParamId{ParameterType::WILSON, "BPrimeCoefficients_B_SCALE_STANDARD", {3051313, 4234, 0, 2}},
        // ParamId{ParameterType::WILSON, "BPrimeCoefficients_B_SCALE_STANDARD", {3051313, 4234, 1, 2}},
        // ParamId{ParameterType::WILSON, "BPrimeCoefficients_B_SCALE_STANDARD", {3051313, 4234, 2, 2}},
        // TODO : Wilsons C10, CQ1, CQ2 + primes at all orders
    }},
    {DecayMapper::to_id(Decays::B__l_nu), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "VCKM", LhaID(0, 2)},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::FLAVOR, "FCONST", LhaID(521,1)}
        // TODO : Wilsons C_V_12, C_S_12 in WGroup::CC_bu
    }},
    {DecayMapper::to_id(Decays::Bs__phi_l_l), {
        
    }},
    {DecayMapper::to_id(Decays::B__Xs_gamma), {
        
    }},
    {DecayMapper::to_id(Decays::B__Xs_l_l), {
        
    }},
    {DecayMapper::to_id(Decays::D__l_nu), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "VCKM", LhaID(1, 0)},
        ParamId{ParameterType::FLAVOR, "FMASS", 411},
        ParamId{ParameterType::FLAVOR, "FLIFE", 411},
        ParamId{ParameterType::FLAVOR, "FCONST", LhaID(411,1)}
        // TODO : Wilsons C_V_12, C_S_12 in WGroup::CC_cd
    }},
    {DecayMapper::to_id(Decays::Ds__l_nu), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "VCKM", LhaID(1, 1)},
        ParamId{ParameterType::FLAVOR, "FMASS", 431},
        ParamId{ParameterType::FLAVOR, "FLIFE", 431},
        ParamId{ParameterType::FLAVOR, "FCONST", LhaID(431,1)}
        // TODO : Wilsons C_V_12, C_S_12 in WGroup::CC_cs
    }},
    {DecayMapper::to_id(Decays::K__l_l), {
        
    }},
    {DecayMapper::to_id(Decays::K__l_nu), {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "VCKM", LhaID(0, 1)},
        ParamId{ParameterType::SM, "VCKM", LhaID(0, 0)},
        ParamId{ParameterType::FLAVOR, "FMASS", 321},
        ParamId{ParameterType::FLAVOR, "FLIFE", 321},
        ParamId{ParameterType::FLAVOR, "FCONST", LhaID(321,1)},
        ParamId{ParameterType::FLAVOR, "FMASS", 211},
        ParamId{ParameterType::FLAVOR, "FLIFE", 211},
        ParamId{ParameterType::FLAVOR, "FCONST", LhaID(211,1)}
        // TODO : Wilsons C_V_12, C_S_12 in WGroup::CC_us and CC_ud
    }},
    {DecayMapper::to_id(Decays::K__pi_nu_nu), {
        
    }},
    {DecayMapper::to_id(Decays::Lambda_b__Lambda_l_l), {
        
    }},
    {DecayMapper::to_id(Decays::M0_Mix), {
        
    }}
};

std::unordered_set<ParamId> DependenciesHelper::get_allowed_parameters(Observables id) {
    ObservableId obs_id = ObservableMapper::to_id(id);

    return get_allowed_parameters(obs_id);
}

bool DependenciesHelper::is_param_allowed(Observables id, ParamId pid) {
    ObservableId obs_id = ObservableMapper::to_id(id);

    return is_param_allowed(obs_id, pid);
}

std::unordered_set<ParamId> DependenciesHelper::get_allowed_parameters(ObservableId id) {
    std::optional<DecayId> did = DecayMapper::get_decay_id(id);
    if (!did.has_value())
        LOG_INFO("No decay found for obs", id.str());
    LOG_INFO("Decay of obs", id.str(), ":", did.value().str());
    if (did.has_value()) {
        return dep_lists.at(did.value());
    } else {
        LOG_ERROR("ValueError", "Observable", ObservableMapper::str(id), "doesn't belong to any declared decay.");
    }
}

bool DependenciesHelper::is_param_allowed(ObservableId id, ParamId pid) {
    auto allowed = get_allowed_parameters(id);
    return std::find(allowed.begin(), allowed.end(), pid) != allowed.end();
}
