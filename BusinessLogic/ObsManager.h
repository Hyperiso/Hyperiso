#ifndef __OBSERVABLEMANAGER_H__
#define __OBSERVABLEMANAGER_H__

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

    double evaluate(Observables id);
    std::map<Observables, double> evaluate_all();
    void add_obs_dep(Observables id, ParamId param);
    void add_obs_deps(Observables id, std::vector<ParamId> params);
    void add_all_obs_deps(Observables id);
    double get_uncertainty(Observables id);
    std::map<Observables, double> get_all_uncertainties();
    std::map<ParamId, double> get_leading_uncertainties(Observables id, size_t n);
    double get_chi2();
    std::vector<Observables> get_current_obss();
    size_t get_obs_evals(Observables id);
    void update_gradient(Observables id);

private:
    ObsManager();
    static std::shared_ptr<ObsManager> instance;

    std::map<Decays, std::shared_ptr<DecayParent>> decays;
    std::map<Observables, std::shared_ptr<Observable>> obss;
    ModelEvaluator me;

    Observables ensure_present(Observables id);

};

#endif // __OBSERVABLEMANAGER_H__
