#include "ObsManager.h"


ObsManager::ObsManager(ObservablePortsConfig obs_port_conf, bool init_default_decays) : obs_port_conf(obs_port_conf) {
    auto& sm = *obs_port_conf.iobspp_sm;
    double mu_W = sm(ParamId{ParameterType::WILSON, "EW_SCALE", 1}, DataType::VALUE);
    double mu_b = sm(ParamId{ParameterType::WILSON, "B_SCALE", 1}, DataType::VALUE);
    if (init_default_decays) {
        this->decays = {
            {DecayMapper::to_id(Decays::B__D_l_nu),             std::make_shared<BDlnuDecay>        (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Dstar_l_nu),         std::make_shared<BDstarlnuDecay>    (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Kstar_gamma),        std::make_shared<BKstarGammaDecay>  (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__l_l),                std::make_shared<BllDecay>          (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__l_nu),               std::make_shared<BlnuDecay>         (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
            {DecayMapper::to_id(Decays::B__Xs_gamma),           std::make_shared<BXsDecay>          (QCDOrder::NONE, mu_W, mu_b, obs_port_conf)},
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

ObsManager ObsManager::add_obs(BinnedObservableId id, QCDOrder order, bool add_deps) {
    const std::string obs_name = ObservableMapper::str(id.s);

    LOG_INFO(
        "Adding binned observable",
        obs_name,
        "with bin [",
        id.p.first,
        ",",
        id.p.second,
        "]"
    );

    auto maybe_obs_enum = ObservableMapper::enum_of(id.s);
    if (!maybe_obs_enum.has_value()) {
        LOG_WARN("Cannot recover enum for observable", obs_name);
        return *this;
    }

    Decays decay_enum;
    try {
        decay_enum = DecayMapper::get_decay(maybe_obs_enum.value());
    } catch (const std::exception& e) {
        LOG_WARN(
            "Cannot find decay for binned observable",
            obs_name,
            ":",
            e.what()
        );
        return *this;
    }

    DecayId dec_id = DecayMapper::to_id(decay_enum);

    auto dec_it = this->decays.find(dec_id);
    if (dec_it == this->decays.end()) {
        LOG_WARN(
            "Decay",
            dec_id.str(),
            "not present in ObsManager for observable",
            obs_name
        );
        return *this;
    }

    dec_it->second->set_order(order);

    if (!this->obss.contains(id.s)) {
        auto obs_ptr = std::make_shared<Observable>(
            id.s,
            dec_it->second,
            obs_port_conf.iobspp_sm
        );

        obss.emplace(id.s, obs_ptr);

        if (add_deps) {
            add_all_obs_deps(id.s);
        }
    }

    dec_it->second->add_bin(id.p);

    LOG_INFO(
        "Successfully added bin [",
        id.p.first,
        ",",
        id.p.second,
        "] for observable",
        obs_name,
        "in decay",
        dec_id.str()
    );

    return *this;
}

// ObsManager ObsManager::add_obs(BinnedObservableId id, QCDOrder order, bool add_deps) {
//     const std::string obs_name = ObservableMapper::str(id.s);
//     LOG_INFO(
//         "Adding binned observable",
//         obs_name,
//         "with bin [",
//         id.p.first,
//         ",",
//         id.p.second,
//         "]"
//     );
//     try {
//         if (!this->obss.contains(id.s)) {
//             add_obs(id.s, order, add_deps);
//         }

//         auto maybe_dec_id = DecayMapper::get_decay_id(id.s);
//         if (!maybe_dec_id.has_value()) {
//             LOG_WARN("No decay id for binned observable", obs_name);
//             return *this;
//         }

//         const DecayId dec_id = maybe_dec_id.value();

//         LOG_INFO(
//             "Resolved binned observable",
//             obs_name,
//             "to decay",
//             dec_id.str()
//         );

//         auto dec_it = this->decays.find(dec_id);
//         if (dec_it == this->decays.end()) {
//             LOG_WARN(
//                 "Decay",
//                 dec_id.str(),
//                 "not present in ObsManager for observable",
//                 obs_name
//             );
//             return *this;
//         }

//         LOG_INFO(
//             "Calling add_bin for",
//             obs_name,
//             "in decay",
//             dec_id.str(),
//             "with bin [",
//             id.p.first,
//             ",",
//             id.p.second,
//             "]"
//         );

//         dec_it->second->add_bin(id.p);

//         LOG_INFO(
//             "Successfully added bin [",
//             id.p.first,
//             ",",
//             id.p.second,
//             "] for observable",
//             obs_name
//         );

//     } catch (const std::out_of_range& e) {
//         LOG_WARN(
//             "Skipping bin [",
//             id.p.first,
//             ",",
//             id.p.second,
//             "] for observable",
//             obs_name,
//             "because std::out_of_range was thrown:",
//             e.what()
//         );
//     } catch (const std::exception& e) {
//         LOG_WARN(
//             "Skipping bin [",
//             id.p.first,
//             ",",
//             id.p.second,
//             "] for observable",
//             obs_name,
//             "because exception was thrown:",
//             e.what()
//         );
//     }
//     // if (!this->obss.contains(id.s))
//     //     add_obs(id.s, order, add_deps);

//     // DecayId dec_id = DecayMapper::get_decay_id(id.s).value();
//     // this->decays.at(dec_id)->add_bin(id.p);
//     return *this;
// }

ObsManager ObsManager::remove_obs(ObservableId id)
{
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

std::map<ObservableId, std::vector<ObservableValue>> ObsManager::evaluate_all() {
    this->enable_obs();
    std::map<ObservableId, std::vector<ObservableValue>> all_ests;
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
    // auto allowed = DependenciesHelper::get_allowed_parameters(id);
    // LOG_DEBUG("Found", allowed.size(), "allowed parameters for observable", ObservableMapper::str(id));
    // obss.at(ensure_present(id))->add_dependences(allowed);
    try {
        auto allowed = DependenciesHelper::get_allowed_parameters(id);
        LOG_DEBUG("Found", allowed.size(), "allowed parameters for observable", ObservableMapper::str(id));
        obss.at(ensure_present(id))->add_dependences(allowed);
    } catch (const std::out_of_range&) {
        LOG_WARN("No dependency map for observable", ObservableMapper::str(id),
                 "- adding it without automatic dependencies.");
    }
}

std::unordered_set<ParamId> ObsManager::get_all_ops_deps(ObservableId id) {
    return DependenciesHelper::get_allowed_parameters(id);
}

std::vector<BinnedObservableId> ObsManager::get_current_obss() {
    std::vector<BinnedObservableId> ids;
    for (const auto& [oid, _] : this->obss) {
        DecayId dec_id = DecayMapper::get_decay_id(oid).value();
        auto bins = this->decays.at(dec_id)->get_bins();
        if (!bins.has_value()) {
            ids.emplace_back(oid);
            continue;
        }
        for (auto bin : bins.value())
            ids.emplace_back(oid, bin);
    }
    return ids;
}

std::shared_ptr<Observable> ObsManager::get_obs(Observables id){
    ObservableId obs_id = ObservableMapper::to_id(id);
    return get_obs(obs_id);
}

std::shared_ptr<Observable> ObsManager::get_obs(ObservableId id){
    return this->obss.at(ensure_present(id));
}

//TODO : better than this
void ObsManager::select_decay(ObservableId id) {
    for (auto& [dec_id, decay] : this->decays) {
        dec_id == DecayMapper::to_id(DecayMapper::get_decay(ObservableMapper::enum_of(id).value())) ? decay->enable() : decay->disable();
    }
}


// void ObsManager::select_decay(ObservableId id) {
//     auto wanted = DecayMapper::to_id(DecayMapper::get_decay(ObservableMapper::enum_of(id).value()));
//     if (active_decay.has_value() && active_decay.value() == wanted) return;

//     if (active_decay.has_value()) {
//         decays.at(active_decay.value())->disable();
//     }

//     decays.at(wanted)->enable();
//     active_decay = wanted;
// }


void ObsManager::set_decay_config(Decays dec, std::any config) {
    this->decays.at(DecayMapper::to_id(dec))->set_config(config);
}

ObservableId ObsManager::ensure_present(Observables id, bool) {
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
    std::unordered_set<DecayId> unique_decays;
    for (auto& elem: this->obss) {
        auto id = DecayMapper::get_decay_id(elem.first);
        if (!id.has_value()) {
            LOG_ERROR("ValueError", "DecayId does not exist for", elem.first.str());
        }
        unique_decays.emplace(id.value());
    }

    for (auto& did : unique_decays)
        this->decays[did]->enable();
}


ObsManager ObsManager::set_bkstarll_threads(size_t n_threads) {
    // LOG_ERROR("NotImplementedError", "No multithreading implemented in old version");
    auto dec_id = DecayMapper::to_id(Decays::B__Kstar_l_l);
    auto dec = std::dynamic_pointer_cast<BKstarllDecay>(this->decays.at(dec_id));
    if (!dec) {
        LOG_ERROR("TypeError", "Decay B__Kstar_l_l is not a BKstarllDecay");
    }
    dec->set_n_threads(n_threads);
    return *this;
}