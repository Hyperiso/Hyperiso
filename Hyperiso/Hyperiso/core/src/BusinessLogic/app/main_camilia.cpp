#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp;
    Config config;
    config.model = Model::CUSTOM;
    config.flags[ExternalFlag::USE_MARTY] = true;
    config.mty_model_name = "ZPrime";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/ZPrime.h");
    hyp.init("lha/camilia.flha", config);
    
    LOG_INFO("HyperisoMaster initialized");

    BlockProxy().log_block(ParameterType::SM, "MASS");
    BlockProxy().log_block(ParameterType::BSM, "MASS");

    auto obs_int = ObservableInterface();

    obs_int.add_observable(Observables::BR_BS_MUMU_UNTAG, QCDOrder::LO, true);
    obs_int.add_observable(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, QCDOrder::LO, true);

    std::cout << obs_int.compute_observable(Observables::BR_BS_MUMU_UNTAG) << std::endl;
    std::cout << obs_int.compute_uncertainty(Observables::BR_BS_MUMU_UNTAG, UncertaintyType::COMBINED) << std::endl;
    std::cout << obs_int.compute_observable(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA) << std::endl;
    std::cout << obs_int.compute_uncertainty(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, UncertaintyType::COMBINED) << std::endl;

    std::cout << obs_int.compute_chi2() << std::endl;

    ParameterSetter().mutate({ParameterType::BSM, "MASS", 32}, 10);

    LOG_INFO(ParameterProvider()({ParameterType::BSM, "MASS", 32}));

    std::cout << obs_int.compute_observable(Observables::BR_BS_MUMU_UNTAG) << std::endl;
    std::cout << obs_int.compute_uncertainty(Observables::BR_BS_MUMU_UNTAG, UncertaintyType::COMBINED) << std::endl;
    std::cout << obs_int.compute_observable(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA) << std::endl;
    std::cout << obs_int.compute_uncertainty(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, UncertaintyType::COMBINED) << std::endl;

    std::cout << obs_int.compute_chi2() << std::endl;
    return 0;
}