#ifndef OBSERVABLEMANAGER_H
#define OBSERVABLEMANAGER_H

#include <memory>
#include <map>

#include "General.h"
#include "Observable.h"
#include "ModelEvaluator.h"
#include "Decays.h"

struct ConfigSetter {
    template <typename EnumType>
    static void apply(DecayParent* base, Decays id, EnumType flag) {
        switch(id) {
            case Decays::B__Kstar_l_l: {
                auto* cfg = static_cast<BKstarllDecay*>(base);
                cfg->set_config_flag(flag);
                break;
            }
            default:
                throw std::logic_error("This decay does not support knobs");
        }
    }
};

class ObsManager {
public:
    ObsManager(std::shared_ptr<ObsWilsonBuilder>& wil_builder);

    ObsManager add_obs(Observables id, QCDOrder order, bool add_deps=false);
    ObsManager remove_obs(Observables id);

    scalar_t evaluate(Observables id);
    std::unordered_map<Observables, Estimate> evaluate_all();
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
    std::shared_ptr<Observable> get_obs(Observables id);

    template <typename EnumType>
    void set_config_flag(Decays decay_id, EnumType flag) {
        if (!decays.contains(decay_id))
            throw std::logic_error("Decay not found in manager");
        ConfigSetter::apply(decays.at(decay_id).get(), decay_id, flag);
    }

    void disable_decays();

private:
    std::unordered_map<Decays, std::shared_ptr<DecayParent>> decays;
    std::unordered_map<Observables, std::shared_ptr<Observable>> obss;
    ModelEvaluator me;

    Observables ensure_present(Observables id, bool critical=true);
    std::shared_ptr<ObsWilsonBuilder> wil_builder;
};

#endif // __OBSERVABLEMANAGER_H__
