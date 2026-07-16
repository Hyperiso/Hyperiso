#include <complex>
#include <iostream>
#include <unordered_set>

#include "HyperisoMaster.h"
#include "ObservableInterface.h"
#include "mapper_hub.hpp"

int main() {

    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);

    GroupMapper::register_custom("DEV_OBS_WILSON", {"dev-obs-wilson"});
    WCoefMapper::register_custom("C_OBS_DEV", {"cobsdev"}, std::pair<int,int>{991001, 1});

    WGroupId group = GroupMapper::id_of("DEV_OBS_WILSON");
    WCoefId coeff = WCoefMapper::id_of("C_OBS_DEV");

    CustomWilsonCoefficientConfig cobs(coeff);
    cobs.set_matching(
        QCDOrder::LO,
        {},
        [](const ParamSrc&) { return 0.75; },
        ContributionType::SM
    );

    CustomWilsonGroupConfig wc(group);
    wc.matching_scale = 160.0;
    wc.hadronic_scale = 4.8;
    wc.order = QCDOrder::LO;
    wc.contribution = ContributionType::SM;
    wc.add_coefficient(cobs);

    LambdaDecayConfig decay;
    decay.canonical = "DEV_LAMBDA_DECAY";
    decay.aliases = {"dev-decay"};
    decay.matching_scale = 160.0;
    decay.hadronic_scale = 4.8;
    decay.order = QCDOrder::LO;
    decay.custom_wilson_groups.push_back(wc);

    LambdaObservableConfig obs = LambdaObservableConfig::scalar(
        "DEV_LAMBDA_OBS",
        [group, coeff](LambdaDecay& ctx, ObservableId id) {
            auto c = ctx.W().getFM(group, coeff, QCDOrder::LO, ContributionType::SM);
            (void)id;
            return std::real(2.0 * c);
        }
    );
    obs.aliases = {"dev-obs"};
    obs.flha = LhaID(992001, 1);
    decay.observables.push_back(obs);

    ObservableInterface O;
    O.add_lambda_decay(decay, true);

    ObservableId obs_id = ObservableMapper::id_of("DEV_LAMBDA_OBS");
    auto values = O.compute_observable(obs_id);

    for (const auto& v : values) {
        std::cout << ObservableMapper::str(v.id) << " = " << v.value << "\n";
    }

    return 0;
}
