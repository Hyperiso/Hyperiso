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
    {DecayMapper::to_id(Decays::B__Xs), {
        
    }},
    {DecayMapper::to_id(Decays::B__Xs_l_l), {
        
    }},
    {DecayMapper::to_id(Decays::D__l_nu), {
        
    }},
    {DecayMapper::to_id(Decays::Ds__l_nu), {
        
    }},
    {DecayMapper::to_id(Decays::K__l_l), {
        
    }},
    {DecayMapper::to_id(Decays::K__l_nu), {
        
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
    if (did.has_value()) {
        dep_lists.at(did.value());
    } else {
        LOG_ERROR("ValueError", "Observable", ObservableMapper::str(id), "doesn't belong to any declared decay.");
    }
}

bool DependenciesHelper::is_param_allowed(ObservableId id, ParamId pid) {
    auto allowed = get_allowed_parameters(id);
    return std::find(allowed.begin(), allowed.end(), pid) != allowed.end();
}
