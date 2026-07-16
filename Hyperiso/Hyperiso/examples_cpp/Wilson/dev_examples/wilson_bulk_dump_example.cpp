#include <iostream>
#include <unordered_set>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "WilsonInterface.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Dev example: dump a full Wilson group at matching and hadronic scales.
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);

    WilsonInterface wilson;
    wilson.build(WilsonBuildConfig(std::unordered_set<WGroup>{WGroup::B}, 81.0, 4.8, QCDOrder::NNLO));

    std::cout << "\nFull matching coefficients, B group, TOTAL, up to NNLO:\n";
    for (const auto& [coef, value] : wilson.getAFM(WGroup::B, QCDOrder::NNLO, ContributionType::TOTAL)) {
        std::cout << "  " << WCoefMapper::str(coef) << " = " << value << "\n";
    }

    std::cout << "\nFull running coefficients, B group, TOTAL, B_STANDARD basis, up to NNLO:\n";
    for (const auto& [coef, value] : wilson.getAFR(WGroup::B, QCDOrder::NNLO, ContributionType::TOTAL, WilsonBasis::B_STANDARD)) {
        std::cout << "  " << WCoefMapper::str(coef) << " = " << value << "\n";
    }

    std::cout << "\nFull running coefficients, B group, TOTAL, B_TRADITIONAL basis, up to NNLO:\n";
    for (const auto& [coef, value] : wilson.getAFR(WGroup::B, QCDOrder::NNLO, ContributionType::TOTAL, WilsonBasis::B_TRADITIONAL)) {
        std::cout << "  " << WCoefMapper::str(coef) << " = " << value << "\n";
    }

    return 0;
}
