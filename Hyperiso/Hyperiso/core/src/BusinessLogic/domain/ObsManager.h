#ifndef OBSERVABLEMANAGER_H
#define OBSERVABLEMANAGER_H

#include <memory>
#include <map>

#include "ParamID.h"
#include "Observable.h"
#include "Decays.h"
#include "IObsParameterProxy.h"
#include "ObsPortsConfig.h"

class ObsManager {
public:
    ObsManager(ObservablePortsConfig obs_port_conf);

    ObsManager add_obs(Observables id, QCDOrder order, bool add_deps=false);
    ObsManager remove_obs(Observables id);

    ObsManager add_obs(ObservableId id, QCDOrder order, bool add_deps=false);
    ObsManager remove_obs(ObservableId id);

    std::vector<ObservableValue> evaluate(Observables id);
    std::vector<ObservableValue> evaluate(ObservableId id);
    std::unordered_map<ObservableId, std::vector<ObservableValue>> evaluate_all();

    void add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr);

    void add_obs_dep(Observables id, ParamId param);
    void add_obs_deps(Observables id, std::unordered_set<ParamId> params);
    void add_obs_dep(ObservableId id, ParamId param);
    void add_obs_deps(ObservableId id, std::unordered_set<ParamId> params);

    std::unordered_set<ParamId> get_all_ops_deps(ObservableId id);

    void add_all_obs_deps(Observables id);

    void add_all_obs_deps(ObservableId id);

    std::unordered_set<ObservableId> get_current_obss();


    std::shared_ptr<Observable> get_obs(Observables id);

    std::shared_ptr<Observable> get_obs(ObservableId id);

    ObservablePortsConfig& get_ports() {return obs_port_conf;}
    
    void select_decay(ObservableId id);

    void reload_params();

    void enable_obs();
    void set_decay_config(Decays dec, std::any config);

private:
    std::unordered_map<DecayId, std::shared_ptr<DecayParent>> decays;
    std::unordered_map<ObservableId, std::shared_ptr<Observable>> obss;
    // ModelEvaluator me;

    ObservableId ensure_present(Observables id, bool critical=true);
    ObservableId ensure_present(ObservableId id, bool critical=true);

    ObservablePortsConfig obs_port_conf;
};

#endif // OBSERVABLEMANAGER_H
