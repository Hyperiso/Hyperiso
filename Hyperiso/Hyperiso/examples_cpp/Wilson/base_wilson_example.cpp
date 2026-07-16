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
    config_hyp.flags[ExternalFlag::HYP_AS_SM_MARTY] = false;
    config_hyp.model = Model::SM;
    config_hyp.mty_model_name = "marty_model_name";
    config_hyp.mty_model_path = "/my/custom/marty/path";

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);

    // Build configuration for Wilson coefficients, with matching/running scales, Wilson groups and order in QCD.
    WilsonBuildConfig config(
        std::unordered_set<WGroup>{WGroup::B, WGroup::BScalar}, // Groups to be built.
        81.0,                         // Matching scale. This scale is global for all groups in this config.
        2.0,                          // Running scale. Special scales for K/D/etc. can be modified with ParameterSetter.
        QCDOrder::LO                  // Order in QCD. With MARTY, SM can go up to NNLO but BSM part stays at LO.
    );

    // Wilson Interface for building and requesting Wilson coefficients.
    WilsonInterface interface;

    // API to build the Wilson coefficients specified within the configuration.
    interface.build(config);

    // You can add a new group after the build. The same kind of QCDOrder and scales are used from this config.
    WilsonBuildConfig config_prime(std::unordered_set<WGroup>{WGroup::BPrime}, 81.0, 2.0, QCDOrder::LO);
    interface.addWilsonGroup(config_prime);

    // Requests to the interface for Wilson coefficients:
    // group: name of the group the coefficient is in.
    // coefficient: name of the coefficient, which needs to be inside the group.
    // order: QCD order at which the coefficient is requested. If above the build order, the coefficient will be 0.
    // contribution: SM, BSM or TOTAL = SM + BSM.
    // scale_type: matching coefficients use getM/getFM; hadronic running coefficients use getR/getFR.
    // WilsonBasis: basis used for the running coefficients. B_STANDARD by default, B_TRADITIONAL also available for B group.
    const std::vector<WCoef> coefs = {
        WCoef::C1, WCoef::C2, WCoef::C3, WCoef::C4, WCoef::C5,
        WCoef::C6, WCoef::C7, WCoef::C8, WCoef::C9, WCoef::C10
    };

    const std::vector<WCoef> coefs_primes = {
        WCoef::CP1, WCoef::CP2, WCoef::CP3, WCoef::CP4, WCoef::CP5, WCoef::CP6,
        WCoef::CP7, WCoef::CP8, WCoef::CP9, WCoef::CP10, WCoef::CPQ1_MU, WCoef::CPQ2_MU
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

    std::cout << "\nAll BScalar matching coefficients at LO:\n";
    for (const auto& [coef, value] : interface.getAM(WGroup::BScalar, QCDOrder::LO, ContributionType::TOTAL)) {
        std::cout << "  " << WCoefMapper::str(coef) << " = " << value << "\n";
    }

    return 0;
}
