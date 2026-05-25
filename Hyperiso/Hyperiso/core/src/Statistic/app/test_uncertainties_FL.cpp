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
#include "ObservableInterfaceProxy.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "StatParamSourcesProxy.h"
#include "StatDependencyPruner.h"
#include "FitAbstraction.h"
#include "NuisanceReader.h"
#include "DefaultNuisancePathsProvider.h"


int main() {
    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;
    hyp.init("lha/si_input.flha", config_hyp);

    auto oint = std::make_shared<ObservableInterface>();

    std::vector<std::pair<double, double>> lhcb_bins = {
        {0.06, 0.98},
        {1.1, 2.5},
        {2.5, 4.0},
        {4.0, 6.0},
        {6.0, 8.0},
        {15.0, 17.0},
        {17.0, 19.0}
    };
    
    // double q2_min = 0.05;
    // double q2_max = 8.20;
    // int N = 500;
    // double bin_width = (q2_max - q2_min) / N;
    // double bin_low = q2_min;

    // for (size_t i = 0; i < N; i++) {
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     bin_low += bin_width; 
    // }

    // q2_min = 14.70;
    // q2_max = 19.00;
    // N = 250;
    // bin_width = (q2_max - q2_min) / N;
    // bin_low = q2_min;

    // for (size_t i = 0; i < N; i++) {
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO);
    //     oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {bin_low, bin_low + bin_width}}, QCDOrder::NNLO); 
    //     bin_low += bin_width; 
    // }

    for (auto& bin: lhcb_bins) {
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO);
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO); 
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_1_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO); 
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO); 
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO); 
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO); 
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO);
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO);
        oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), bin}, QCDOrder::NNLO);
    }

    BKstarllConfig cfg;
    cfg.ff_src = BV_FF_Src::GRvDV;
    oint->set_decay_config(Decays::B__Kstar_l_l, cfg);
    oint->set_bkstarll_threads(25);
    std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();

    auto model = std::make_shared<ObservableInterfaceProxy>(oint, spop);

    StatisticConfig config;
    config.MC_draws = 1000;
    config.nuisance_relevance_cutoff = 1e-10;
    config.override_nuisance_marginals = {
        {ParamId{ParameterType::WILSON, "EW_SCALE", 1}, MarginalType::FLAT},
        {ParamId{ParameterType::WILSON, "B_SCALE", 1}, MarginalType::FLAT},
    };

    std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();

    StatisticManager stat(
        config,
        model,
        std::make_shared<StatCorrelationProxy>(),
        std::make_shared<StatParameterProxy>(),
        std::make_shared<StatParamSourcesProxy>(),
        std::make_shared<StatDependencyPruner>(),
        std::make_shared<NuisanceReader>(npp),
        spop
    );

    // stat.select_experiment("CMS");

    auto start = std::chrono::steady_clock::now();
    auto pred_with_u = stat.compute_uncertainties();
    auto stop = std::chrono::steady_clock::now();
    LOG_INFO("Uncertainty computation took", std::chrono::duration_cast<std::chrono::seconds>(stop - start).count(), "s");

    std::ofstream fs;
    std::map<std::string, Observables> file_names = {
        {"dGamma_dq2", Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU},
        {"F_L", Observables::F_L_B0__KSTAR0_MU_MU},
        {"P_1", Observables::P_1_B0__KSTAR0_MU_MU},
        {"P_2", Observables::P_2_B0__KSTAR0_MU_MU},
        {"P_3", Observables::P_3_B0__KSTAR0_MU_MU},
        {"Pp_4", Observables::P_PRIME_4_B0__KSTAR0_MU_MU},
        {"Pp_5", Observables::P_PRIME_5_B0__KSTAR0_MU_MU},
        {"Pp_6", Observables::P_PRIME_6_B0__KSTAR0_MU_MU},
        {"Pp_8", Observables::P_PRIME_8_B0__KSTAR0_MU_MU}
    };

    for (auto &&[name, obs] : file_names) {
        fs.open(name + "_LHCB_bins.csv");
        fs << "bin_low,bin_high,value,u_low,u_high,u_sym\n";

        for (auto &&[k, v] : pred_with_u) {
            if (k.s != ObservableMapper::to_id(obs)) continue;
            fs << k.p.first << "," << k.p.second << "," << v.mu << "," << v.sigma_m << "," << v.sigma_p << "," << v.sigma << '\n';
        }
        fs.close();
    }

    return 0;
}