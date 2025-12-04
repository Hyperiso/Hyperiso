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
    }

    std::map<ObservableId, double> compute_uncertainties() {
        return manager->compute_uncertainties();
    }
    
    FitResultWithMaps compute_MLE() {
        return manager->compute_MLE();
    }

private:
    std::shared_ptr<StatisticManager> manager;
};