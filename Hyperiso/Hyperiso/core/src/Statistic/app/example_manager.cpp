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

    config.p_specs = {ParamId(ParameterType::SM, "MASS", 1)};
    std::shared_ptr<ObservableInterface> oi = std::make_shared<ObservableInterface>();

    StatisticManager stat(config, std::make_shared<ObservableInterfaceAdapterObs>(oi), std::make_shared<StatCorrelationProxy>(), std::make_shared<StatParameterProxy>(), std::make_shared<StatParamSourcesProxy>());

    stat.fill_cache();

    auto start = std::chrono::steady_clock::now();
    stat.compute_uncertainties();
    auto stop  = std::chrono::steady_clock::now();

    auto us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    std::cout << "Temps : " << us << " µs\n";

    auto res = stat.compute_MLE();

    std::cout << "Temps : " << res.ell_hat << " µs\n";

    ParamId pid = config.p_specs[0];

    // auto interval = stat.compute_CI_1d_95(pid, -7.0, -1.0, 50);

    auto interval = stat.compute_CI_1d_95(pid, 0.001, 0.09, 10);

    std::cout << "IC 95% sur " << pid << " = ["
            << interval.first << ", " << interval.second << "]\n";

    // std::cout << "Temps : " << res.eta_hat[0] << " µs\n";
    // std::cout << "Temps : " << res.p_hat[0] << " µs\n";

    return 0;
}
