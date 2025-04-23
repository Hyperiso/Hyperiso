#include "ObsManager.h"

std::shared_ptr<ObsManager> ObsManager::instance = nullptr;

ObsManager::ObsManager() {
    ObsParameterProxy obsParamProxy = ObsParameterProxy(ParameterType::FLAVOR);
    

    this->decays = {
        {Decays::B__D_l_nu, std::make_shared<BDlnuDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 521))},
        {Decays::B__Dstar_l_nu, std::make_shared<BDstarlnuDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 521))},
        {Decays::B__Kstar,  std::make_shared<BKstarDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 511))},
        {Decays::B__l_l,    std::make_shared<BllDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 531))},
        {Decays::B__l_nu,   std::make_shared<BlnuDecay>(QCDOrder::NONE, 81, obsParamProxy("FMASS", 511))},
        {Decays::B__Xs,     std::make_shared<BXsDecay>(QCDOrder::NONE, 81, ObsParameterProxy(ParameterType::SM)("QCD", LhaID(5,3)) / 2)},
    };

    // this->decays = {
    //     {Decays::B__D_l_nu, std::make_shared<BDlnuDecay>(QCDOrder::NONE, 81, Parameters::Get(ParameterType::FLAVOR, "FMASS", 521))},
    //     {Decays::B__Dstar_l_nu, std::make_shared<BDstarlnuDecay>(QCDOrder::NONE, 81, Parameters::Get(ParameterType::FLAVOR, "FMASS", 521))},
    //     {Decays::B__Kstar,  std::make_shared<BKstarDecay>(QCDOrder::NONE, 81, Parameters::Get(ParameterType::FLAVOR, "FMASS", 511))},
    //     {Decays::B__l_l,    std::make_shared<BllDecay>(QCDOrder::NONE, 81, Parameters::Get(ParameterType::FLAVOR, "FMASS", 531))},
    //     {Decays::B__l_nu,   std::make_shared<BlnuDecay>(QCDOrder::NONE, 81, Parameters::Get(ParameterType::FLAVOR, "FMASS", 511))},
    //     {Decays::B__Xs,     std::make_shared<BXsDecay>(QCDOrder::NONE, 81, QCDHelper::mass_b_1S() / 2)},
    // };
}

std::shared_ptr<ObsManager> ObsManager::GetInstance() {
    if (!instance) {
        instance = std::shared_ptr<ObsManager>(new ObsManager());
    }
    return instance;
}

std::shared_ptr<ObsManager> ObsManager::add_obs(Observables id, QCDOrder order, bool add_deps) {
    auto dec = decays.at(DecayMapper::get_decay(id));
    dec->set_order(order);
    auto obs_ptr = std::make_shared<Observable>(id, dec);
    obss.emplace(id, obs_ptr);
    me.add_observable(obs_ptr);

    if (add_deps) {
        add_all_obs_deps(id);
    }

    return instance;
}

std::shared_ptr<ObsManager> ObsManager::remove_obs(Observables id) {
    id = ensure_present(id);
    obss.erase(id);
    me.remove_observable(id);

    return instance;
}

double ObsManager::evaluate(Observables id) {
    return obss.at(ensure_present(id))->eval();
}

std::map<Observables, double> ObsManager::evaluate_all() {
    std::map<Observables, double> all_vals;
    for (auto &[k, _] : obss) {
        all_vals.emplace(k, evaluate(k));
    }
    return all_vals;
}

void ObsManager::add_obs_dep(Observables id, ParamId param) {
    if (DependenciesHelper::is_param_allowed(id, param)) {
        obss.at(ensure_present(id))->add_dependence(param);
    }
    else {
        LOG_WARN("Observable", ObservableMapper::str(id), "doesn't depend on parameter (", param.block, ",", param.code, "). Ignoring.");
    }
}

void ObsManager::add_obs_deps(Observables id, std::vector<ParamId> params) {
    for (auto &p : params) {
        if (!DependenciesHelper::is_param_allowed(id, p)) {
            LOG_WARN("Observable", ObservableMapper::str(id), "doesn't depend on parameter (", p.block, ",", p.code, "). Ignoring.");
        }
    }
    obss.at(ensure_present(id))->add_dependences(params);
}

void ObsManager::add_all_obs_deps(Observables id) {
    obss.at(ensure_present(id))->add_dependences(DependenciesHelper::get_allowed_parameters(id));
}

double ObsManager::get_uncertainty(Observables id) {
    return std::sqrt(obss.at(ensure_present(id))->variance());
}

std::map<Observables, double> ObsManager::get_all_uncertainties() {
    std::map<Observables, double> all_vars;
    for (auto &[k, _] : obss) {
        all_vars.emplace(k, get_uncertainty(k));
    }
    return all_vars;
}

std::map<ParamId, double> ObsManager::get_leading_uncertainties(Observables id, size_t n) {
    return obss.at(ensure_present(id))->get_leading_uncertainties(n);
}

double ObsManager::get_chi2() {
    return me.chi2();
}

size_t ObsManager::get_obs_evals(Observables id) {
    return obss.at(ensure_present(id))->get_n_evals();
}

void ObsManager::update_gradient(Observables id) {
    obss.at(ensure_present(id))->update_gradient();
}

Observables ObsManager::ensure_present(Observables id) {
    if (!obss.contains(id)) {
        LOG_ERROR("KeyError", "Observable manager doesn't contain observable", ObservableMapper::str(id));
    }
    return id;
}
