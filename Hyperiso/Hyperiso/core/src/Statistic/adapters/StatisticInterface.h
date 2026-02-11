#pragma once
#include "StatisticManager.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"

class StatisticInterface {
public:
    StatisticInterface(StatisticConfig config) {
        std::shared_ptr<ObservableInterface> oi = std::make_shared<ObservableInterface>();
        std::shared_ptr<IModel> oia= std::make_shared<ObservableInterfaceAdapterObs>(oi);
        std::shared_ptr<IStatCorrelationProxy> pscp = std::make_shared<StatCorrelationProxy>();
        std::shared_ptr<IStatParameterProxy> pspp = std::make_shared<StatParameterProxy>();
        std::shared_ptr<IStatSourcesProxy> sp = std::make_shared<StatParamSourcesProxy>();
        manager = std::make_shared<StatisticManager>(config, oia, pscp, pspp, sp);
        manager->fill_cache();
    }

    std::map<ObservableId, GaussianSummary> compute_uncertainties() {
        return manager->compute_uncertainties();
    }

    MCResult compute_uncertainties_and_sampling() {
        return manager->compute_uncertainties_and_sampling();
    }
    
    // FitResultWithMaps compute_MLE() {
    //     return manager->compute_MLE();
    // }

private:
    std::shared_ptr<StatisticManager> manager;
};