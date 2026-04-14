#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
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

    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;
    hyp.init("lha/si_input.flha", config_hyp);

    auto oint = std::make_shared<ObservableInterface>();

    oint->add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
        .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true);

    auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint);

    StatisticConfig config;
    config.MLE_max_iter = 120000;
    config.MLE_tol = 0.2;

    std::shared_ptr<INuisancePathsProvider> npp =
        std::make_shared<DefaultNuisancePathsProvider>();

    StatisticManager stat(
        config,
        model,
        std::make_shared<StatCorrelationProxy>(),
        std::make_shared<StatParameterProxy>(),
        std::make_shared<StatParamSourcesProxy>(),
        std::make_shared<StatDependencyPruner>(),
        std::make_shared<NuisanceReader>(npp)
    );

    std::vector<ParamId> p_specs = {
        ParamId{ParameterType::WILSON,
                GroupMapper::str(WGroup::B, ScaleType::HADRONIC),
                WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::TOTAL)},
        ParamId{ParameterType::WILSON,
                GroupMapper::str(WGroup::B, ScaleType::HADRONIC),
                WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::TOTAL)}
    };

    stat.prepare_likelihood_for_scan(p_specs);

    std::map<ParamId, double> p_hat_manual = {
        {p_specs[0], 3.8586816584445423},
        {p_specs[1], -4.342464564024878}
    };


    std::map<ParamId, double> eta_hat_manual = {
        {ParamId{ParameterType::WILSON, "B_SCALE", LhaID(1)}, 7.0616412623634615},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,1)}, -1.3194750934885218},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,2)},  0.55731952538812446},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,1)},  0.42794295488555201},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,2)},  1.3824396326751143},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,1)},  0.55297187610318355},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,2)},  0.63064666791673996},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,1)}, -1.1020600964364089},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,2)},  2.0650881994252286},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,1)}, -0.97034002393979002},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,2)},  1.3955789220459496},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,1)},  0.52904322583131957},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,2)},  1.6397824945976849},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,1)},  1.3404825108445939},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,2)},  3.4992796210570716},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(7,1)}, 0.046984916688010275},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(7,2)}, 0.10307659688624525},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(8,1)}, 0.060283567240217828},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(8,2)}, 0.16002917156326671},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(14)},  0.47361707939582515},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,1)},  0.096779421539748903},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,2)},  0.03548473538478969},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,3)},  0.022403500826101146},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,4)},  0.026932380671486327},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,5)}, -0.011150073481602584},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,6)},  0.016123830166622971},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,1)},  0.36090946885144704},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,2)},  0.13795059765754347},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,3)},  0.16977339600345981},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,4)},  0.11216066278714854},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,5)}, -0.060833674931900818},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,6)},  0.046546022227100338},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,1)}, -0.018041968289152464},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,2)},  0.0008849972829968639},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,3)}, -0.063630276500820357},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,4)},  0.00018648240700748546},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,5)},  0.062387459629891652},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,6)},  0.00062480775226799617},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,7)},  5.2618056451280428e-07},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,8)}, -0.099938779409649459}
    };


    stat.set_manual_scan_point(p_hat_manual, eta_hat_manual);

    auto grid = stat.scan_likelihood_around_current_point(
        p_specs[0],
        p_specs[1],
        2.0,
        2.5,
        81,
        81
    );

    stat.save_likelihood_scan_csv("scan_from_manual_point.csv", grid);
    std::cout << "[INFO] Wrote scan_from_manual_point.csv\n";

    return 0;
}
