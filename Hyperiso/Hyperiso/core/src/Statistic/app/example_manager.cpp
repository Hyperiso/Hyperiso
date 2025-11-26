#include "StatisticManager.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"

int main(int argc, char** argv) {

    HyperisoMaster hyp = HyperisoMaster();
    Config config_hyp;
    config_hyp.model = Model::SM;

    hyp.init("lha/si_input.flha", config_hyp);

    StatisticConfig config;

    config.obss = {
        {ObservableMapper::to_id(Observables::BR_BS_MUMU), QCDOrder::LO},
        {ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO},
        {ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO}
    };

    std::shared_ptr<ObservableInterface> oi = std::make_shared<ObservableInterface>();

    StatisticManager stat(config, std::make_shared<ObservableInterfaceAdapterObs>(oi), std::make_shared<StatCorrelationProxy>(), std::make_shared<StatParameterProxy>(), std::make_shared<StatParamSourcesProxy>());

    stat.fill_cache();

    return 0;
}
