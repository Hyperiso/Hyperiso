#pragma once
#include "StatisticManager.h"
#include "ObservableInterfaceAdapter2.h"
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
        std::shared_ptr<IModel> oia= std::make_shared<ObservableInterfaceAdapterObs>(oi);
        std::shared_ptr<IStatCorrelationProxy> pscp = std::make_shared<StatCorrelationProxy>();
        std::shared_ptr<IStatParameterProxy> pspp = std::make_shared<StatParameterProxy>();
        std::shared_ptr<IStatSourcesProxy> sp = std::make_shared<StatParamSourcesProxy>();
        std::shared_ptr<IStatDependencyPruner> sdp = std::make_shared<StatDependencyPruner>();
        std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();
        std::shared_ptr<INuisanceReader> nr = std::make_shared<NuisanceReader>(npp);
        manager = std::make_shared<StatisticManager>(config, oia, pscp, pspp, sp, sdp, nr);
        manager->update_cache();
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