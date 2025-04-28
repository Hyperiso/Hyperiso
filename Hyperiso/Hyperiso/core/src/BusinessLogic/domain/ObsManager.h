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
    static std::shared_ptr<ObsManager> GetInstance();

    std::shared_ptr<ObsManager> add_obs(Observables id, QCDOrder order, bool add_deps=false);
    std::shared_ptr<ObsManager> remove_obs(Observables id);

    scalar_t evaluate(Observables id);
    std::unordered_map<Observables, scalar_t> evaluate_all();
    void add_obs_dep(Observables id, ParamId param);
    void add_obs_deps(Observables id, std::unordered_set<ParamId> params);
    void add_all_obs_deps(Observables id);
    scalar_t get_uncertainty(Observables id);
    std::unordered_map<Observables, scalar_t> get_all_uncertainties();
    std::unordered_map<ParamId, scalar_t> get_leading_uncertainties(Observables id, size_t n);
    double get_chi2();
    std::unordered_set<Observables> get_current_obss();
    size_t get_obs_evals(Observables id);
    void update_gradient(Observables id);

private:
    ObsManager();
    static std::shared_ptr<ObsManager> instance;

    std::unordered_map<Decays, std::shared_ptr<DecayParent>> decays;
    std::unordered_map<Observables, std::shared_ptr<Observable>> obss;
    ModelEvaluator me;

    Observables ensure_present(Observables id);

};

#endif // __OBSERVABLEMANAGER_H__
