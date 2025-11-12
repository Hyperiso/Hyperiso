#include "StatCorrelationProxy.h"
#include "HyperisoMaster.h"
#include "BlockProxy.h"
#include "StatParameterProxy.h"

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

    std::cout << "Corr between A_3_B__KSTAR_L_L and Observables::A_4_B__KSTAR_L_L : "<< stat_corr(Observables::A_3_B__KSTAR_L_L, Observables::A_4_B__KSTAR_L_L, StatCorrelationProxy::Type::COMBINED) << std::endl;

    std::cout << "Corr between same param : "<< stat_corr(Observables::A_3_B__KSTAR_L_L, Observables::A_3_B__KSTAR_L_L, StatCorrelationProxy::Type::COMBINED) << std::endl;

    BlockProxy().log_all_blocks(ParameterType::OBSERVABLE);

    StatParameterProxy spp = StatParameterProxy(ParameterType::SM);

    std::cout << "BR_BD_MUMU obs param : " << *spp.get_obs_param(ObservableMapper::to_id(Observables::BR_BD_MUMU)) << std::endl;

    std::cout << "W mass : " << *spp.get_param(ParamId(ParameterType::SM, "MASS", 24)) << std::endl;
    return 0;
}
