#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "StatisticManager.h"
#include "ObservableInterface.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "StatParamSourcesProxy.h"
#include "StatDependencyPruner.h"
#include "FitAbstraction.h"
#include "NuisanceReader.h"
#include "DefaultNuisancePathsProvider.h"


int main() {
    using namespace fit_app;
    // Logger::getInstance()->setLevel(Logger::LogLevel::VERBOSE);
    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;
    hyp.init("lha/si_input.flha", config_hyp);

    auto oint = std::make_shared<ObservableInterface>();

    double q2_min = 0.05;
    double q2_max = 8.00;
    int N = 100;
    double bin_width = (q2_max - q2_min) / N;
    double bin_low = q2_min;

    for (size_t i = 0; i < N; i++) {
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);  
        bin_low += bin_width; 
    }

    BKstarllConfig cfg;
    cfg.ff_src = BV_FF_Src::GRvDV;
    oint->set_decay_config(Decays::B__Kstar_l_l, cfg);

    auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint);

    StatisticConfig config;
    config.MC_draws = 1000;

    LOG_INFO("Creating StatisticManager.");

    std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();

    StatisticManager stat(
        config,
        model,
        std::make_shared<StatCorrelationProxy>(),
        std::make_shared<StatParameterProxy>(),
        std::make_shared<StatParamSourcesProxy>(),
        std::make_shared<StatDependencyPruner>(),
        std::make_shared<NuisanceReader>(npp)
    );

    LOG_INFO("StatisticManager created.");

    auto start = std::chrono::steady_clock::now();
    auto pred_with_u = stat.compute_uncertainties();
    auto stop = std::chrono::steady_clock::now();

    std::ofstream fs;
    fs.open("F_L.csv");
    fs << "bin_low,bin_high,value,u_low,u_high,u_sym\n";

    for (auto &&[k, v] : pred_with_u) {
        fs << k.p.first << "," << k.p.second << "," << v.mu << "," << v.sigma_m << "," << v.sigma_p << "," << v.sigma << '\n';
    }

    LOG_INFO("Uncertainty computation took", std::chrono::duration_cast<std::chrono::seconds>(stop - start).count(), "s");

    return 0;
}