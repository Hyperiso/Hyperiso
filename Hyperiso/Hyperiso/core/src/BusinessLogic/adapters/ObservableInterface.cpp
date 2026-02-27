#include "ObservableInterface.h"

ObservableInterface::ObservableInterface() {
    WilsonBuildConfig config;
    std::shared_ptr<WilsonBuilder> builder_ptr = std::make_shared<WilsonBuilder>(config);
    std::shared_ptr<IObsWilsonBuilder> builder = std::make_shared<ObsWilsonBuilder>(builder_ptr);
    

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm = std::make_shared<ObsParameterProxy>();

    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav = std::make_shared<ObsParameterProxy>(ParameterType::FLAVOR);
    std::shared_ptr<IObsQCDProxy> iobs_qcdp = std::make_shared<ObsQCDProxy>();


    std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty = std::make_shared<ObsUseMarty>();

    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer = std::make_shared<WilsonFreezer>(builder);
    
    ports = std::make_shared<ObservablePortsConfig>(builder, iobspp_sm, iobspp_flav, iobs_qcdp, iobs_use_marty, iobs_wfreezer);

    manager = std::make_shared<ObsManager>(*ports);
}

void ObservableInterface::add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr) {
    manager->add_custom_decay(id, ptr);
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

void ObservableInterface::add_observables(Decays decay, QCDOrder order, bool add_dependencies) {  
    for (auto &obs : DecayMapper::get_observables(decay)) {
        add_observable(obs, order, add_dependencies);
    }
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
    return manager->get_obs(id)->get_exp_val();
};

scalar_t ObservableInterface::get_exp_value(ObservableId id) {
    return manager->get_obs(id)->get_exp_val();
};

scalar_t ObservableInterface::get_exp_uncertainty(Observables id, UncertaintyType u_type) {
    return manager->get_obs(id)->get_exp_uncertainty(u_type);
}

scalar_t ObservableInterface::get_exp_uncertainty(ObservableId id, UncertaintyType u_type) {
    return manager->get_obs(id)->get_exp_uncertainty(u_type);
}

std::vector<BinnedObservableId> ObservableInterface::get_current_observables() {
    return manager->get_current_obss();
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