#ifndef OBSERVABLEMANAGER_H
#define OBSERVABLEMANAGER_H

#include <memory>
#include <map>

#include "General.h"
#include "Observable.h"
#include "ModelEvaluator.h"
#include "Decays.h"

class ObsManager {
public:
    ObsManager(std::shared_ptr<ObsWilsonBuilder>& wil_builder);

    ObsManager add_obs(Observables id, QCDOrder order, bool add_deps=false);
    ObsManager remove_obs(Observables id);

    ObsManager add_obs(ObservableId id, QCDOrder order, bool add_deps=false);
    ObsManager remove_obs(ObservableId id);

    std::vector<ObservableValue> evaluate(Observables id);
    std::vector<ObservableValue> evaluate(ObservableId id);
    std::unordered_map<ObservableId, Estimate> evaluate_all();

    void add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr);

    void add_obs_dep(Observables id, ParamId param);
    void add_obs_deps(Observables id, std::unordered_set<ParamId> params);
    void add_obs_dep(ObservableId id, ParamId param);
    void add_obs_deps(ObservableId id, std::unordered_set<ParamId> params);

    std::unordered_set<ParamId> get_all_ops_deps(ObservableId id);

    void add_all_obs_deps(Observables id);
    scalar_t get_uncertainty(Observables id);

    void add_all_obs_deps(ObservableId id);
    scalar_t get_uncertainty(ObservableId id);

    std::unordered_map<ObservableId, scalar_t> get_all_uncertainties();
    std::unordered_map<ParamId, scalar_t> get_leading_uncertainties(Observables id, size_t n);

    std::unordered_map<ParamId, scalar_t> get_leading_uncertainties(ObservableId id, size_t n);

    double get_chi2();
    std::unordered_set<ObservableId> get_current_obss();
    void update_gradient(Observables id);

    void update_gradient(ObservableId id);

    std::shared_ptr<Observable> get_obs(Observables id);

    std::shared_ptr<Observable> get_obs(ObservableId id);

    void select_decay(ObservableId id);

private:
    std::unordered_map<DecayId, std::shared_ptr<DecayParent>> decays;
    std::unordered_map<ObservableId, std::shared_ptr<Observable>> obss;
    ModelEvaluator me;

    ObservableId ensure_present(Observables id, bool critical=true);
    ObservableId ensure_present(ObservableId id, bool critical=true);

    std::shared_ptr<ObsWilsonBuilder> wil_builder;
};

#endif // __OBSERVABLEMANAGER_H__
