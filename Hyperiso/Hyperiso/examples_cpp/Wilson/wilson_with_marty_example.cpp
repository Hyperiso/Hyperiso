#include <iostream>
#include <vector>
#include <unordered_set>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "WilsonInterface.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // See base_core_example.cpp for informations about these classes.
    HyperisoConfig config_hyp;
    config_hyp.flags[ExternalFlag::IS_LHA_SPECTRUM] = false;
    config_hyp.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
    config_hyp.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT] = false;
    config_hyp.flags[ExternalFlag::HYP_AS_SM_MARTY] = true;

    // We specify here that we want to use MARTY to calculate the Wilson coefficients.
    config_hyp.model = Model::MARTY;
    config_hyp.mty_model_name = "SM";
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/sm.h";

    HyperisoMaster hyp;

    // Be sure to have MARTY installed. If you installed it within Hyperiso, you don't have anything to do.
    // If you have MARTY on your machine please use the following line:
    // hyp.pre_init_set_marty_path("/path/to/marty");

    // If you have in your LHA new blocks which are not in the FLHA convention, please add this line:
    // hyp.pre_init_add_block("NAME_OF_THE_BLOCK");
    // By default blocks are id + value. If you have another type of block, use the options of pre_init_add_block.

    hyp.init("lha/zprime_input.flha", config_hyp);

    // The rest of the code is the same than the base_wilson_example.cpp example.
    WilsonBuildConfig config(std::unordered_set<WGroup>{WGroup::B}, 81.0, 2.0, QCDOrder::LO);

    WilsonInterface interface;
    interface.build(config);

    WilsonBuildConfig config_prime(std::unordered_set<WGroup>{WGroup::BPrime}, 81.0, 2.0, QCDOrder::LO);
    interface.addWilsonGroup(config_prime);

    const std::vector<WCoef> coefs = {
        WCoef::C1, WCoef::C2, WCoef::C3, WCoef::C4, WCoef::C5,
        WCoef::C6, WCoef::C7, WCoef::C8, WCoef::C9, WCoef::C10
    };

    const std::vector<WCoef> coefs_primes = {
        WCoef::CP1, WCoef::CP2, WCoef::CP3, WCoef::CP4, WCoef::CP5, WCoef::CP6,
        WCoef::CP7, WCoef::CP8, WCoef::CP9, WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2
    };

    std::cout << "\nSeparated matching orders for B group:\n";
    for (auto coef : coefs) {
        std::cout << WCoefMapper::str(coef) << "\n";
        for (const auto& [order, value] : interface.getSM(WGroup::B, coef, ContributionType::TOTAL)) {
            std::cout << "  " << OrderMapper::str(order) << " = " << value << "\n";
        }
    }

    std::cout << "\nFull running coefficients for BPrime group:\n";
    for (auto coef : coefs_primes) {
        std::cout << WCoefMapper::str(coef) << " = "
                  << interface.getFR(WGroup::BPrime, coef, QCDOrder::LO, ContributionType::TOTAL) << "\n";
    }

    return 0;
}
