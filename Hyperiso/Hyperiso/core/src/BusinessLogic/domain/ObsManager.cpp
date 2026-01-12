#include "ObsManager.h"


ObsManager::ObsManager(ObservablePortsConfig obs_port_conf, bool init_default_decays) : obs_port_conf(obs_port_conf) {


    auto& sm = *obs_port_conf.iobspp_sm;
    auto& flav = *obs_port_conf.iobspp_flav;

    double mu_W = 2. * sm(ParamId{ParameterType::SM, "MASS", 24}, DataType::VALUE); // 2 * m_W
    double mu_b = sm(ParamId{ParameterType::SM, "QCD", LhaID(5, 2)}, DataType::VALUE) / 2; // m_b(pole) / 2
    if (init_default_decays) {
        this->decays = {
            {DecayMapper::to_id(Decays::B__D_l_nu),             std::make_shared<BDlnuDecay>        (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Dstar_l_nu),         std::make_shared<BDstarlnuDecay>    (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Kstar_gamma),        std::make_shared<BKstarGammaDecay>  (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__l_l),                std::make_shared<BllDecay>          (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__l_nu),               std::make_shared<BlnuDecay>         (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Xs),                 std::make_shared<BXsDecay>          (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Xs_l_l),             std::make_shared<BXsllDecay>        (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::M0_Mix),                std::make_shared<M0Mixing>          (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Kstar_l_l),          std::make_shared<BKstarllDecay>     (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__K_l_l),              std::make_shared<BKllDecay>         (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::Bs__phi_l_l),           std::make_shared<BsPhiDecay>        (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::Lambda_b__Lambda_l_l),  std::make_shared<LbLllDecay>        (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::K__l_l),                std::make_shared<KllDecay>          (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::K__pi_nu_nu),           std::make_shared<KPinunuDecay>      (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::K__l_nu),               std::make_shared<KlnuDecay>         (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::D__l_nu),               std::make_shared<DlnuDecay>         (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::Ds__l_nu),              std::make_shared<DslnuDecay>        (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
        };
    }
}

ObsManager ObsManager::add_obs(ObservableId id, QCDOrder order, bool add_deps) {
    LOG_INFO("Adding observable", ObservableMapper::str(id), "to manager");
    // LOG_INFO(DecayMapper::get_decay_id(id).value().str()); //TODO : error with get_decay_id
    auto dec = decays.at(DecayMapper::to_id(DecayMapper::get_decay(ObservableMapper::enum_of(id).value())));
    dec->set_order(order);
    auto obs_ptr = std::make_shared<Observable>(id, dec, obs_port_conf.iobspp_sm);
    obss.emplace(id, obs_ptr);
    // me.add_observable(obs_ptr);
    if (add_deps) {
        add_all_obs_deps(id);
    }
    return *this;
}

ObsManager ObsManager::add_obs(Observables id, QCDOrder order, bool add_deps) {
    ObservableId obs_id = ObservableMapper::to_id(id);
    return add_obs(obs_id, order, add_deps);
}

ObsManager ObsManager::remove_obs(Observables id) {
    ObservableId obs_id = ObservableMapper::to_id(id);
    return remove_obs(obs_id);
}

ObsManager ObsManager::remove_obs(ObservableId id) {
    id = ensure_present(id, false);
    obss.erase(id);
    // me.remove_observable(id);

    return *this;
}

std::vector<ObservableValue> ObsManager::evaluate(ObservableId id) {
    select_decay(id);
    return obss.at(ensure_present(id))->compute();
}

std::vector<ObservableValue> ObsManager::evaluate(Observables id) {
    ObservableId obs_id = ObservableMapper::to_id(id);
    return evaluate(obs_id);
}

std::unordered_map<ObservableId, std::vector<ObservableValue>> ObsManager::evaluate_all() {
    std::unordered_map<ObservableId, std::vector<ObservableValue>> all_ests;
    for (auto &[k, k_ptr] : obss) {
        all_ests.emplace(k, k_ptr->compute());
    }
    return all_ests;
}

void ObsManager::add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr) {
    ptr->bind_wilson_builder(this->obs_port_conf.iobswb);
    // ptr->load_params();
    this->decays.emplace(id, ptr);
}

void ObsManager::add_obs_dep(Observables id, ParamId param) {
    ObservableId obs_id = ObservableMapper::to_id(id);

    return add_obs_dep(obs_id, param);
}

void ObsManager::add_obs_dep(ObservableId id, ParamId param) {
    if (DependenciesHelper::is_param_allowed(id, param)) {
        obss.at(ensure_present(id))->add_dependence(param);
    }
    else {
        LOG_WARN("Observable", ObservableMapper::str(id), "doesn't depend on parameter (", param.block, ",", param.code, "). Ignoring.");
    }
}

void ObsManager::add_obs_deps(Observables id, std::unordered_set<ParamId> params) {
    ObservableId obs_id = ObservableMapper::to_id(id);

    return add_obs_deps(obs_id, params);
}

void ObsManager::add_obs_deps(ObservableId id, std::unordered_set<ParamId> params) {
    for (auto &p : params) {
        if (ObservableMapper::enum_of(id).has_value() && !DependenciesHelper::is_param_allowed(id, p)) {
            LOG_WARN("Observable", ObservableMapper::str(id), "doesn't depend on parameter (", p.block, ",", p.code, "). Ignoring.");
        }
    }
    obss.at(ensure_present(id))->add_dependences(params);
}

void ObsManager::add_all_obs_deps(Observables id) {
    ObservableId obs_id = ObservableMapper::to_id(id);

    return add_all_obs_deps(obs_id);
}

void ObsManager::add_all_obs_deps(ObservableId id) {
    auto allowed = DependenciesHelper::get_allowed_parameters(id);
    LOG_DEBUG("Found", allowed.size(), "allowed parameters for observable", ObservableMapper::str(id));
    obss.at(ensure_present(id))->add_dependences(allowed);
}

std::unordered_set<ParamId> ObsManager::get_all_ops_deps(ObservableId id) {
    return DependenciesHelper::get_allowed_parameters(id);
}

std::unordered_set<ObservableId> ObsManager::get_current_obss() {
    return get_keys(this->obss);
}

std::shared_ptr<Observable> ObsManager::get_obs(Observables id){
    ObservableId obs_id = ObservableMapper::to_id(id);
    return get_obs(obs_id);
}

std::shared_ptr<Observable> ObsManager::get_obs(ObservableId id){
    return this->obss.at(ensure_present(id));
}

void ObsManager::select_decay(ObservableId id) {
    for (auto& [dec_id, decay] : this->decays) {
        dec_id == DecayMapper::to_id(DecayMapper::get_decay(ObservableMapper::enum_of(id).value())) ? decay->enable() : decay->disable();
    }
}

void ObsManager::set_decay_config(Decays dec, std::any config) {
    this->decays.at(DecayMapper::to_id(dec))->set_config(config);
}

ObservableId ObsManager::ensure_present(Observables id, bool critical) {
    ObservableId obs_id = ObservableMapper::to_id(id);
    return ensure_present(obs_id);
}

ObservableId ObsManager::ensure_present(ObservableId id, bool critical) {
    if (!obss.contains(id)) {
        if (critical)
            LOG_ERROR("KeyError", "Observable manager doesn't contain observable", ObservableMapper::str(id));
        else
            LOG_WARN("Observable manager doesn't contain observable", ObservableMapper::str(id));
    }
    return id;
}

void ObsManager::reload_params() {
    for (auto& elem: this->decays) {
        elem.second->load_params();
    }
}

void ObsManager::enable_obs() {
    for (auto& elem: this->obss) {
        auto id = DecayMapper::get_decay_id(elem.first);
        if (!id.has_value()) {
            LOG_ERROR("ValueError", "DecayId does not exist for", elem.first.str());
        }
        this->decays[id.value()]->disable();
        this->decays[id.value()]->enable(); //TODO : avoid double enabling (two obs of same decay) please
    }
}