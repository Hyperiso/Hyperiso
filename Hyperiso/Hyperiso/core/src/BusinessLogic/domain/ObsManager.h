#ifndef OBSERVABLEMANAGER_H
#define OBSERVABLEMANAGER_H

#include <memory>
#include <map>

#include "ParamID.h"
#include "Observable.h"
#include "Decays.h"
#include "IObsParameterProxy.h"

struct ObservablePortsConfig {
    /**
     * @brief Constructs the port bundle.
     *
     * @param iblock_c         Block composer used to register dependent blocks/parameters.
     * @param wilson_proxy     Proxy to read/write Wilson parameters by (block, LhaID).
     * @param use_marty        Runtime flag: whether Marty backend is enabled.
     * @param has_wilson       Runtime flag: whether an input FWCOEF block is present.
     * @param model_api        Runtime access to the currently selected physics model.
     * @param scale_setter_api Runtime setter for switching and setting the active scale.
     */
    ObservablePortsConfig(std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> iobswb, std::shared_ptr<IObsWilsonProxy<ObsWilsonBuilder>> iobswp, std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp, std::shared_ptr<IObsQCDProxy> iobs_qcdp, std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty, std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer) :
        iobswb(iobswb), iobswp(iobswp), 
        iobspp(iobspp), iobs_qcdp(iobs_qcdp),
        iobs_use_marty(iobs_use_marty), iobs_wfreezer(iobs_wfreezer) {}

    /// Dependency engine used to compose derived blocks/parameters.
    std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> iobswb;

    /// Read-only access to existing blocks/parameters (FWCOEF, matching blocks, etc.).
    std::shared_ptr<IObsWilsonProxy<ObsWilsonBuilder>> iobswp;

    /// Backend selector: true => Marty, false => builtin.
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp;

    /// Whether an input Wilson block exists (FWCOEF-style input present).
    std::shared_ptr<IObsQCDProxy> iobs_qcdp;

    /// Current model selector (SM/SUSY/THDM...).
    std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty;

    /// Scale setter used to switch and set mu_W / mu_h.
    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer;

    /**
     * @brief Optional group builder hook.
     *
     * Signature: (group id, model, use_marty, contribution type, storage block name) -> group instance.
     *
     * When provided, the manager uses it to build:
     *  - SM-only intermediate groups (stored in a block name like "<MATCHING>_SM"),
     *  - possibly other specializations depending on the application.
     *
     * If null, the manager falls back to cloning existing groups via @ref CoefficientGroup::get_sm_group().
     */
    std::function<std::shared_ptr<CoefficientGroup>(WGroupId, Model, bool, ContributionType, std::string)> build_group;
};

class ObsManager {
public:
    ObsManager(std::shared_ptr<ObsWilsonBuilder>& wil_builder);

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

    std::shared_ptr<ObsWilsonBuilder> wil_builder;
};

#endif // OBSERVABLEMANAGER_H
