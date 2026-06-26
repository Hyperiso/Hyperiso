#include "ObservableInterface.h"

ObservableInterface::ObservableInterface() {
    WilsonBuildConfig config;
    std::shared_ptr<WilsonBuilder> builder_ptr = std::make_shared<WilsonBuilder>(config);
    std::shared_ptr<IObsWilsonBuilder> builder = std::make_shared<ObsWilsonBuilder>(builder_ptr);
    std::shared_ptr<ObsWilsonHelper> wilson_helper = std::make_shared<ObsWilsonHelper>();
    

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm = std::make_shared<ObsParameterProxy>();

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav = std::make_shared<ObsParameterProxy>(ParameterType::FLAVOR);
    std::shared_ptr<IObsQCDProxy> iobs_qcdp = std::make_shared<ObsQCDProxy>();


    std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty = std::make_shared<ObsUseMarty>();

    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer = std::make_shared<WilsonFreezer>(builder);
    
    ports = std::make_shared<ObservablePortsConfig>(builder, iobspp_sm, iobspp_flav, iobs_qcdp, iobs_use_marty, iobs_wfreezer, wilson_helper);

    manager = std::make_shared<ObsManager>(*ports);
}

void ObservableInterface::add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr) {
    manager->add_custom_decay(id, ptr);
}

ObservableInterface& ObservableInterface::add_lambda_decay(LambdaDecayConfig config, bool add_observables) {
    if (config.observables.empty()) {
        LOG_ERROR("ValueError", "Cannot register lambda decay", config.canonical, "without observables.");
    }

    std::vector<CustomObservableSpec> specs;
    specs.reserve(config.observables.size());
    for (const auto& obs : config.observables) {
        specs.push_back(CustomObservableSpec{obs.canonical, obs.aliases, obs.flha});
    }

    const bool registered = DecayMapper::register_custom_with_observables(
        config.canonical,
        config.aliases,
        std::move(specs)
    );

    if (!registered) {
        LOG_WARN(
            "Lambda decay",
            config.canonical,
            "was not newly registered. Reusing an existing mapper entry if it exists."
        );
    }

    for (const auto& custom_group : config.custom_wilson_groups) {
        this->ports->iobswb->add_custom_group(custom_group);
        this->ports->iobs_whelper->mark_built(
            custom_group.group,
            custom_group.matching_scale,
            custom_group.hadronic_scale,
            custom_group.order
        );
    }

    DecayId decay_id = DecayMapper::id_of(config.canonical);
    auto decay = std::make_shared<LambdaDecay>(decay_id, config, *this->ports);
    this->manager->add_custom_decay(decay_id, decay);

    if (add_observables) {
        for (const auto& obs : config.observables) {
            ObservableId obs_id = ObservableMapper::id_of(obs.canonical);
            this->manager->add_obs(obs_id, config.order, false);
            if (!obs.dependencies.empty()) {
                this->manager->add_obs_deps(obs_id, obs.dependencies);
            }
        }
    }

    return *this;
}

ObservableInterface& ObservableInterface::add_observable(Observables obs, QCDOrder order, bool add_dependencies) {
    manager->add_obs(obs, order, add_dependencies);
    return *this;
}

ObservableInterface& ObservableInterface::add_observable(ObservableId obs, QCDOrder order, bool add_dependencies) { 
    manager->add_obs(obs, order, add_dependencies);
    return *this;
}

ObservableInterface &ObservableInterface::add_observable(BinnedObservableId obs,
                                                         QCDOrder order,
                                                         bool add_dependencies)
{
    manager->add_obs(obs, order, add_dependencies);
    return *this;
}

void ObservableInterface::add_observables(std::map<Observables, QCDOrder> obss,
                                          bool add_dependencies)
{
    for (auto &[k, v] : obss) {
        add_observable(k, v, add_dependencies);
    }
}

void ObservableInterface::add_observables(std::map<ObservableId, QCDOrder> obss, bool add_dependencies) {  
    for (auto &[k, v] : obss) {

        add_observable(k, v, add_dependencies);
    }
}

void ObservableInterface::add_observables(Decays decay,
                                          QCDOrder order,
                                          bool add_dependencies,
                                          std::pair<double, double> bin)
{
    for (auto &obs : DecayMapper::get_observables(decay)) {
        if (manager->is_observable_binned(obs)) {
            add_observable(BinnedObservableId(obs, bin), order, add_dependencies);
        } else {
            add_observable(obs, order, add_dependencies);
        }
    }
}

void ObservableInterface::add_observables(DecayId decay,
                                          QCDOrder order,
                                          bool add_dependencies,
                                          std::pair<double, double> bin)
{
    for (auto &obs : DecayMapper::get_observables(decay)) {
        if (manager->is_observable_binned(obs)) {
            add_observable(BinnedObservableId(obs, bin), order, add_dependencies);
        } else {
            add_observable(obs, order, add_dependencies);
        }
    }
}

bool ObservableInterface::is_decay_binned(Decays decay) const {
    return manager->is_decay_binned(decay);
}

bool ObservableInterface::is_decay_binned(DecayId decay) const {
    return manager->is_decay_binned(decay);
}

bool ObservableInterface::is_observable_binned(Observables obs) const {
    return manager->is_observable_binned(obs);
}

bool ObservableInterface::is_observable_binned(ObservableId obs) const {
    return manager->is_observable_binned(obs);
}

void ObservableInterface::add_observable_parameter(Observables obs, ParamId pid) {
    manager->add_obs_dep(obs, pid);
}

void ObservableInterface::add_observable_parameter(ObservableId obs, ParamId pid) {
    manager->add_obs_dep(obs, pid);
}

void ObservableInterface::add_observable_parameters(Observables obs, std::unordered_set<ParamId> pids) {
    manager->add_obs_deps(obs, pids);
}

void ObservableInterface::add_observable_parameters(ObservableId obs, std::unordered_set<ParamId> pids) {
    manager->add_obs_deps(obs, pids);
}

std::vector<ObservableValue> ObservableInterface::compute_observable(Observables obs) const {
    return manager->evaluate(obs);
}

std::vector<ObservableValue> ObservableInterface::compute_observable(ObservableId obs) const {
    return manager->evaluate(obs);
}

ObservableValue ObservableInterface::compute_observable(BinnedObservableId obs) const {
    return manager->evaluate(obs);
}

void ObservableInterface::remove_observable(Observables id) {
    manager->remove_obs(id);
}

void ObservableInterface::remove_observable(ObservableId id) {
    manager->remove_obs(id);
}

void ObservableInterface::remove_observables(std::unordered_set<Observables> ids) {
    for (Observables id : ids) {
        remove_observable(id);
    }
};


void ObservableInterface::remove_observables(Decays dec) {
    for (Observables id : DecayMapper::get_observables(dec)) {
        remove_observable(id);
    }
};

scalar_t ObservableInterface::get_exp_value(Observables id) {
    return manager->get_obs(id)->get_exp_val({0,0}, "DEFAULT");
};

scalar_t ObservableInterface::get_exp_value(ObservableId id) {
    return manager->get_obs(id)->get_exp_val({0,0}, "DEFAULT");
};

scalar_t ObservableInterface::get_exp_value(BinnedObservableId id) {
    return manager->get_obs(id.s)->get_exp_val(id.p, "DEFAULT");
};

scalar_t ObservableInterface::get_exp_value(ExperimentObs id) {
    return manager->get_obs(id.obs.s)->get_exp_val(id.obs.p, id.experiment);
};

scalar_t ObservableInterface::get_exp_uncertainty(Observables id, UncertaintyType u_type) {
    return manager->get_obs(id)->get_exp_uncertainty({0,0}, "DEFAULT", u_type);
}

scalar_t ObservableInterface::get_exp_uncertainty(ObservableId id, UncertaintyType u_type) {
    return manager->get_obs(id)->get_exp_uncertainty({0,0}, "DEFAULT", u_type);
}

scalar_t ObservableInterface::get_exp_uncertainty(BinnedObservableId id, UncertaintyType u_type) {
    return manager->get_obs(id.s)->get_exp_uncertainty(id.p, "DEFAULT", u_type);
}

scalar_t ObservableInterface::get_exp_uncertainty(ExperimentObs id, UncertaintyType u_type) {
    return manager->get_obs(id.obs.s)->get_exp_uncertainty(id.obs.p, id.experiment, u_type);
}

std::vector<BinnedObservableId> ObservableInterface::get_current_observables() {
    return manager->get_current_obss();
}

std::shared_ptr<ObservableInterface> ObservableInterface::clone_for_worker() const {
    auto worker = std::make_shared<ObservableInterface>();
    auto selection = manager->snapshot_observable_selection();

    for (const auto& item : selection) {
        if (item.binned) {
            worker->add_observable(item.id, item.order, false);
        } else {
            worker->add_observable(item.id.s, item.order, false);
        }

        if (!item.dependencies.empty()) {
            worker->add_observable_parameters(item.id.s, item.dependencies);
        }
    }

    return worker;
}

DecayThreadSnapshot ObservableInterface::snapshot_decay_threads() const {
    return manager->snapshot_decay_threads();
}

void ObservableInterface::set_all_decay_threads(size_t n_threads) {
    manager->set_all_decay_threads(n_threads);
}

void ObservableInterface::restore_decay_threads(const DecayThreadSnapshot& snapshot) {
    manager->restore_decay_threads(snapshot);
}

std::map<ObservableId, std::vector<ObservableValue>> ObservableInterface::compute_all() {
    return manager->evaluate_all();
}


void ObservableInterface::set_param(const std::string& block, LhaID code, double value, ParameterType type) {
    Parameters::GetInstance(type)->setBlockValue(block, code, value);
}

double ObservableInterface::get_param(const std::string& block, LhaID code, ParameterType type) {
    return Parameters::GetInstance(type)->operator()(block, code);
}

std::unordered_set<ParamId> ObservableInterface::get_all_ops_deps(ObservableId id) {
    return manager->get_all_ops_deps(id);
}

std::unordered_set<ParamId> ObservableInterface::get_all_ops_deps(Observables id) {
    return manager->get_all_ops_deps(ObservableMapper::to_id(id));
}

void ObservableInterface::reload_params() {
    manager->reload_params();
}

void ObservableInterface::enable_obs() {
    manager->enable_obs();
}
void ObservableInterface::set_decay_config(Decays dec, std::any config) {
    manager->set_decay_config(dec, config);  
}

ObservablePortsConfig& ObservableInterface::get_ports() {return this->manager->get_ports();}

void ObservableInterface::set_decay_threads(Decays dec, size_t n_threads) {this->manager->set_decay_threads(dec, n_threads);}

void ObservableInterface::set_bkstarll_threads(size_t n_threads) {this->manager->set_bkstarll_threads(n_threads);}

void ObservableInterface::set_bkll_threads(size_t n_threads) {this->manager->set_bkll_threads(n_threads);}

void ObservableInterface::set_bsphi_threads(size_t n_threads) {this->manager->set_bsphi_threads(n_threads);}

void ObservableInterface::set_lblll_threads(size_t n_threads) {this->manager->set_lblll_threads(n_threads);}