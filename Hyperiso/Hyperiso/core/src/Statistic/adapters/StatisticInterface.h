#pragma once
#include "StatisticManager.h"
#include "ObservableInterfaceProxy.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"
#include "StatDependencyPruner.h"
#include "DefaultNuisancePathsProvider.h"
#include "NuisanceReader.h"

class StatisticInterface {
public:
    StatisticInterface(StatisticConfig config, std::shared_ptr<ObservableInterface> oi_) {
        std::shared_ptr<ObservableInterface> oi = oi_;
        std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
        std::shared_ptr<IModel> oia = std::make_shared<ObservableInterfaceProxy>(oi, spop);
        std::shared_ptr<IStatCorrelationProxy> pscp = std::make_shared<StatCorrelationProxy>();
        std::shared_ptr<IStatParameterProxy> pspp = std::make_shared<StatParameterProxy>();
        std::shared_ptr<IStatSourcesProxy> sp = std::make_shared<StatParamSourcesProxy>();
        std::shared_ptr<IStatDependencyPruner> sdp = std::make_shared<StatDependencyPruner>();
        std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();
        std::shared_ptr<INuisanceReader> nr = std::make_shared<NuisanceReader>(npp);

        manager = std::make_shared<StatisticManager>(config, oia, pscp, pspp, sp, sdp, nr, spop);

        manager->select_experiments_all();
        manager->update_cache();
    }

    void select_experiment(const std::string& experiment) {
        manager->select_experiment(experiment);
        manager->update_cache();
    }

    void select_experiments(const std::vector<std::string>& experiments) {
        manager->select_experiments(experiments);
        manager->update_cache();
    }

    void select_experiments_all() {
        manager->select_experiments_all();
        manager->update_cache();
    }

    bool has_experiment_selection() const noexcept {
        return manager->has_experiment_selection();
    }

    std::set<std::string> selected_experiments() const {
        return manager->selected_experiments();
    }

    std::map<BinnedObservableId, GaussianSummary> compute_uncertainties() {
        return manager->compute_uncertainties();
    }

    MCResult compute_uncertainties_and_sampling() {
        return manager->compute_uncertainties_and_sampling();
    }
    
    FitResultWithMaps compute_MLE(const std::vector<ParamId>& p_specs) {
        return manager->compute_MLE(p_specs);
    }

    Contour compute_confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
        return manager->confidence_contour(p1, p2, z, bounds, options);
    }

private:
    std::shared_ptr<StatisticManager> manager;
};