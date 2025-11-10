#include "StatCorrelationProxy.h"
#include "HyperisoMaster.h"
int main(int argc, char** argv) {

    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::SM;

    hyp.init("lha/testInput.flha", config);

    auto stat_corr = StatCorrelationProxy();

    double el = stat_corr(ObservableMapper::from_flha(LhaID("531_15_2_13_-13")).value(),
    ObservableMapper::from_flha(LhaID("511_1_2_13_-13")).value(),
    StatCorrelationProxy::Type::COMBINED);

    std::cout << el << std::endl;

    std::cout << stat_corr(Observables::A_3_B__KSTAR_L_L, Observables::A_4_B__KSTAR_L_L, StatCorrelationProxy::Type::COMBINED) << std::endl;

    std::cout << stat_corr(Observables::A_3_B__KSTAR_L_L, Observables::A_3_B__KSTAR_L_L, StatCorrelationProxy::Type::COMBINED) << std::endl;
    return 0;
}
