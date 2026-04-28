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

    std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
    auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint, spop);

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
        std::make_shared<NuisanceReader>(npp),
        spop
    );

    const std::string had_bsm_block =
        GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD)
        + "__BSM_INTERMEDIATE";

    std::vector<ParamId> p_specs = {
        ParamId{ParameterType::WILSON,
                had_bsm_block,
                WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, ContributionType::BSM)},
        ParamId{ParameterType::WILSON,
                had_bsm_block,
                WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, ContributionType::BSM)}
    };

    stat.prepare_likelihood_for_scan(p_specs);

    // std::map<ParamId, double> p_hat_manual = {
    //     {p_specs[0], 3.8586816584445423},
    //     {p_specs[1], -4.342464564024878}
    // };
    // std::map<ParamId, double> p_hat_manual = {
    //     {p_specs[0], -0.55},
    //     {p_specs[1], 0}
    // };

    std::map<ParamId, double> p_hat_manual = {
        {p_specs[0], -4.1913054387387545},
        {p_specs[1],  4.0360988362819041}
    };

    // std::map<ParamId, double> eta_hat_manual = {
    //     {ParamId{ParameterType::WILSON, "B_SCALE", LhaID(1)}, 7.0616412623634615},

    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,1)}, -1.3194750934885218},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,2)},  0.55731952538812446},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,1)},  0.42794295488555201},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,2)},  1.3824396326751143},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,1)},  0.55297187610318355},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,2)},  0.63064666791673996},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,1)}, -1.1020600964364089},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,2)},  2.0650881994252286},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,1)}, -0.97034002393979002},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,2)},  1.3955789220459496},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,1)},  0.52904322583131957},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,2)},  1.6397824945976849},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,1)},  1.3404825108445939},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,2)},  3.4992796210570716},

    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(7,1)}, 0.046984916688010275},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(7,2)}, 0.10307659688624525},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(8,1)}, 0.060283567240217828},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(8,2)}, 0.16002917156326671},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(14)},  0.47361707939582515},

    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,1)},  0.096779421539748903},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,2)},  0.03548473538478969},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,3)},  0.022403500826101146},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,4)},  0.026932380671486327},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,5)}, -0.011150073481602584},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,6)},  0.016123830166622971},

    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,1)},  0.36090946885144704},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,2)},  0.13795059765754347},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,3)},  0.16977339600345981},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,4)},  0.11216066278714854},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,5)}, -0.060833674931900818},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,6)},  0.046546022227100338},

    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,1)}, -0.018041968289152464},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,2)},  0.0008849972829968639},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,3)}, -0.063630276500820357},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,4)},  0.00018648240700748546},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,5)},  0.062387459629891652},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,6)},  0.00062480775226799617},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,7)},  5.2618056451280428e-07},
    //     {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,8)}, -0.099938779409649459},

    //     // Nouvelles nuisances Wilson détectées après le fix des blocs intermédiaires
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,4422,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,4422,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,4422,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,6421,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,6421,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,6421,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,4141,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,4141,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,4141,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,6161,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,6161,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,6161,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4133,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4133,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4133,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4536,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4536,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4536,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6153,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6153,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6153,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6556,0,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6556,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6556,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4133,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4133,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4137,1,1)}, 0.0},
    //     {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4137,2,1)}, 0.0},

    //     {ParamId{ParameterType::WILSON, "EW_SCALE", LhaID(1)}, 100.0},
    // };

    std::map<ParamId, double> eta_hat_manual = {
        // MASS
        {ParamId{ParameterType::SM, "MASS", LhaID(1)}, 0.0047000246138147843},
        {ParamId{ParameterType::SM, "MASS", LhaID(2)}, 0.0021600098832398637},
        {ParamId{ParameterType::SM, "MASS", LhaID(3)}, 0.093493634506413584},
        {ParamId{ParameterType::SM, "MASS", LhaID(4)}, 1.2729652178959461},
        {ParamId{ParameterType::SM, "MASS", LhaID(24)}, 80.369199884059867},
        {ParamId{ParameterType::SM, "MASS", LhaID(25)}, 125.20003079073854},

        // SMINPUTS
        {ParamId{ParameterType::SM, "SMINPUTS", LhaID(3)}, 0.11807867791958603},
        {ParamId{ParameterType::SM, "SMINPUTS", LhaID(5)}, 4.1852875835851924},
        {ParamId{ParameterType::SM, "SMINPUTS", LhaID(6)}, 172.56568043942119},
        {ParamId{ParameterType::SM, "SMINPUTS", LhaID(7,1)}, 0.23116000044751273},

        // VCKMIN
        {ParamId{ParameterType::SM, "VCKMIN", LhaID(1)}, 0.22500998835442002},
        {ParamId{ParameterType::SM, "VCKMIN", LhaID(2)}, 0.82599937438110005},
        {ParamId{ParameterType::SM, "VCKMIN", LhaID(3)}, 0.15913636223243458},
        {ParamId{ParameterType::SM, "VCKMIN", LhaID(4)}, 0.35228967135353356},

        // FLIFE
        {ParamId{ParameterType::FLAVOR, "FLIFE", LhaID(511)}, 1.5169999568101041e-12},
        {ParamId{ParameterType::FLAVOR, "FLIFE", LhaID(521)}, 1.6379999957628896e-12},

        // FMASS
        {ParamId{ParameterType::FLAVOR, "FMASS", LhaID(313)}, 0.8955496332249121},

        // Wilson nuisances
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,4422,0,1)}, 0.35408395412360966},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,4422,1,1)}, 0.0051048552231385165},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,4422,2,1)}, 7.260983968910622e-05},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,6421,0,1)}, 0.64251236865612282},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,6421,1,1)}, 0.010054467169024599},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(305,6421,2,1)}, 0.00016288296132363063},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,4141,0,1)}, -0.017214870710872351},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,4141,1,1)}, -0.00081968012872567751},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,4141,2,1)}, -1.3103861690999478e-05},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,6161,0,1)}, 0.73174221890563995},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,6161,1,1)}, 0.011305085833610801},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3040405,6161,2,1)}, 0.00018337129877690338},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4133,0,1)}, -0.030577049689067357},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4133,1,1)}, -0.0022774584267058516},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4133,2,1)}, -7.5077683845400526e-05},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4536,0,1)}, -0.016478446922171562},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4536,1,1)}, -0.0019444051939759893},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,4536,2,1)}, -0.00068495697454101247},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6153,0,1)}, -0.44588339089840923},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6153,1,1)}, -0.0073645682464370957},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6153,2,1)}, -0.00011941893854722827},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6556,0,1)}, 0.059466834352638671},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6556,1,1)}, 0.00050744582582298143},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3050707,6556,2,1)}, 8.7390117885074317e-06},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4133,1,1)}, 0.00021701970443030607},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4133,2,1)}, 3.8650299507191088e-06},

        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4137,1,1)}, 0.00041002202481149543},
        {ParamId{ParameterType::WILSON, had_bsm_block, LhaID(3051313,4137,2,1)}, 7.0954953085568712e-06},

        // scales
        {ParamId{ParameterType::WILSON, "B_SCALE", LhaID(1)}, 5.8851667928051219},
        {ParamId{ParameterType::WILSON, "EW_SCALE", LhaID(1)}, 100.72686727075933},

        // B_Ks
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,0)}, 0.36935806519425035},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,1)}, -1.3769142196051429},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,1,2)}, 0.034633126485268971},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,0)}, 0.29149202558029325},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,1)}, 0.37537257434702137},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,2,2)}, 1.2401493646876649},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,0)}, 0.26552633766960121},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,1)}, 0.53379909757705002},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,3,2)}, 0.45740855009424058},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,0)}, 0.37157858125584586},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,1)}, -1.1935759144599001},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,4,2)}, 2.6121881087516545},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,0)}, 0.30694856227428019},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,1)}, -1.0293275682602907},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,5,2)}, 1.7209514285050549},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,0)}, 0.30693011922964947},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,1)}, 0.48101238727118456},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,6,2)}, 1.6942364262993246},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,0)}, 0.66514149858850979},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,1)}, 1.3070135823523328},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(1,7,2)}, 3.8510126859669498},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(7,1)}, 0.041277395162247661},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(7,2)}, 0.10624518694150206},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(8,1)}, 0.059290556077814471},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(8,2)}, 0.16044517867574637},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,1)}, -0.069778583601311284},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,2)}, -0.0053711735591936479},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,3)}, 0.03647857273215125},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,4)}, 0.039693408296069148},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,5)}, -0.0040592226041542973},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,1,6)}, -0.0063410854580069082},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,1)}, -0.08192548685195418},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,2)}, 0.0057286354696556444},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,3)}, 0.22607364917793951},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,4)}, -0.027886074215018514},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,5)}, -0.0074003272932216639},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,2,6)}, -0.017999786361442286},

        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,1)}, 0.031697199815891496},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,2)}, -0.00027091449357067497},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,3)}, -0.052885117244455694},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,4)}, -0.0052498994138951158},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,5)}, 0.019257463859406802},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,6)}, 0.0032318869717003164},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,7)}, 3.383515277131575e-05},
        {ParamId{ParameterType::DECAY, "B_Ks", LhaID(18,3,8)}, 1.7810438682197634e-05},
    };

    stat.set_manual_scan_point(p_hat_manual, eta_hat_manual);

    auto grid = stat.scan_likelihood_around_current_point(
        p_specs[0],
        p_specs[1],
        10.0,
        10,
        40,
        40
    );

    stat.save_likelihood_scan_csv("scan_from_manual_point.csv", grid);
    std::cout << "[INFO] Wrote scan_from_manual_point.csv\n";

    return 0;
}
