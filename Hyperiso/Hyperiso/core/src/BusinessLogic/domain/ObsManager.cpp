#include "ObsManager.h"


ObsManager::ObsManager(std::shared_ptr<ObsWilsonBuilder>& wil_builder) {
    this->wil_builder = wil_builder;
    ObsParameterProxy obsParamProxy = ObsParameterProxy(ParameterType::FLAVOR);
    ObsParameterProxy smParamProxy = ObsParameterProxy(ParameterType::SM);
    
    this->decays = {
        {Decays::B__D_l_nu,     std::make_shared<BDlnuDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 521), wil_builder)},
        {Decays::B__Dstar_l_nu, std::make_shared<BDstarlnuDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 521), wil_builder)},
        {Decays::B__Kstar,      std::make_shared<BKstarDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 511), wil_builder)},
        {Decays::B__l_l,        std::make_shared<BllDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 531), wil_builder)},
        {Decays::B__l_nu,       std::make_shared<BlnuDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 511), wil_builder)},
        {Decays::B__Xs,         std::make_shared<BXsDecay>(QCDOrder::NONE, 81, smParamProxy("QCD", LhaID(5, 3)) / 2, wil_builder)},
        {Decays::B__Xs_l_l,     std::make_shared<BXsllDecay>(QCDOrder::NONE, 2 * smParamProxy("MASS", 24), smParamProxy("QCD", LhaID(5, 2)) / 2, wil_builder)},
        {Decays::M0_Mix,        std::make_shared<M0Mixing>(QCDOrder::NONE, 160, 4.16, wil_builder)},
        {Decays::B__Kstar_l_l,  std::make_shared<BKstarllDecay>(QCDOrder::NONE, 81, smParamProxy("QCD", LhaID(5, 3)) / 2, wil_builder)},
    };
}

ObsManager ObsManager::add_obs(Observables id, QCDOrder order, bool add_deps) {
    LOG_INFO("Adding observable", ObservableMapper::str(id), "to manager");
    auto dec = decays.at(DecayMapper::get_decay(id));
    dec->set_order(order);
    auto obs_ptr = std::make_shared<Observable>(id, dec);
    obss.emplace(id, obs_ptr);
    me.add_observable(obs_ptr);
    if (add_deps) {
        add_all_obs_deps(id);
    }
    return *this;
}

ObsManager ObsManager::remove_obs(Observables id) {
    id = ensure_present(id, false);
    obss.erase(id);
    me.remove_observable(id);

    return *this;
}

scalar_t ObsManager::evaluate(Observables id) {
    return obss.at(ensure_present(id))->eval();
}

std::unordered_map<Observables, Estimate> ObsManager::evaluate_all() {
    std::unordered_map<Observables, Estimate> all_ests;
    for (auto &[k, k_ptr] : obss) {
        all_ests.emplace(k, k_ptr->get_estimate());
    }
    return all_ests;
}

void ObsManager::add_obs_dep(Observables id, ParamId param) {
    if (DependenciesHelper::is_param_allowed(id, param)) {
        obss.at(ensure_present(id))->add_dependence(param);
    }
    else {
        LOG_WARN("Observable", ObservableMapper::str(id), "doesn't depend on parameter (", param.block, ",", param.code, "). Ignoring.");
    }
}

void ObsManager::add_obs_deps(Observables id, std::unordered_set<ParamId> params) {
    for (auto &p : params) {
        if (!DependenciesHelper::is_param_allowed(id, p)) {
            LOG_WARN("Observable", ObservableMapper::str(id), "doesn't depend on parameter (", p.block, ",", p.code, "). Ignoring.");
        }
    }
    obss.at(ensure_present(id))->add_dependences(params);
}

void ObsManager::add_all_obs_deps(Observables id) {
    auto allowed = DependenciesHelper::get_allowed_parameters(id);
    LOG_DEBUG("Found", allowed.size(), "allowed parameters for observable", ObservableMapper::str(id));
    obss.at(ensure_present(id))->add_dependences(allowed);
}

scalar_t ObsManager::get_uncertainty(Observables id) {
    return std::sqrt(obss.at(ensure_present(id))->variance());
}

std::unordered_map<Observables, scalar_t> ObsManager::get_all_uncertainties() {
    std::unordered_map<Observables, scalar_t> all_vars;
    for (auto &[k, _] : obss) {
        all_vars.emplace(k, get_uncertainty(k));
    }
    return all_vars;
}

std::unordered_map<ParamId, scalar_t> ObsManager::get_leading_uncertainties(Observables id, size_t n) {
    return obss.at(ensure_present(id))->get_leading_uncertainties(n);
}

double ObsManager::get_chi2() {
    return me.chi2();
}

std::unordered_set<Observables> ObsManager::get_current_obss() {
    return get_keys(this->obss);
}

size_t ObsManager::get_obs_evals(Observables id) {
    return obss.at(ensure_present(id))->get_n_evals();
}

void ObsManager::update_gradient(Observables id) {
    obss.at(ensure_present(id))->update_gradient();
}

std::shared_ptr<Observable> ObsManager::get_obs(Observables id){
    return this->obss.at(ensure_present(id));
}

void ObsManager::disable_decays() {
    for (auto& [_, decay] : this->decays) {
        decay->disable();
    }
}

Observables ObsManager::ensure_present(Observables id, bool critical) {
    if (!obss.contains(id)) {
        if (critical)
            LOG_ERROR("KeyError", "Observable manager doesn't contain observable", ObservableMapper::str(id));
        else
            LOG_WARN("Observable manager doesn't contain observable", ObservableMapper::str(id));
    }
    return id;
}
