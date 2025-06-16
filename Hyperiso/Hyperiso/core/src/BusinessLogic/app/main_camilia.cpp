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
    // ParameterSetter().mutate({ParameterType::BSM, "MASS", 32}, 1.5);
    
    LOG_INFO("HyperisoMaster initialized");

    auto obs_int = ObservableInterface();

    obs_int.add_observable(Observables::BR_B_XS_GAMMA, QCDOrder::LO, true);


    std::cout << obs_int.compute_observable(Observables::BR_B_XS_GAMMA) << std::endl;
    std::cout << obs_int.compute_uncertainty(Observables::BR_B_XS_GAMMA, UncertaintyType::COMBINED) << std::endl;

    obs_int.add_observable(Observables::BR_BS_MUMU, QCDOrder::LO);

    std::cout << obs_int.get_exp_value(Observables::BR_BS_MUMU) << " +- " << obs_int.get_exp_uncertainty(Observables::BR_BS_MUMU) << std::endl; 
    obs_int.add_observable(Observables::BR_BD_MUMU, QCDOrder::LO, true);

    std::cout << obs_int.compute_observable(Observables::BR_BD_MUMU) << std::endl;
    std::cout << obs_int.compute_uncertainty(Observables::BR_BD_MUMU, UncertaintyType::COMBINED) << std::endl;

    std::cout << obs_int.compute_chi2() << std::endl;
    return 0;
}