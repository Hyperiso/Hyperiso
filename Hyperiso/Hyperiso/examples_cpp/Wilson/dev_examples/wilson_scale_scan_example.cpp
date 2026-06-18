#include <iostream>
#include <vector>
#include <unordered_set>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "WilsonInterface.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Dev example: scan the hadronic scale without rebuilding the whole interface.
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);

    WilsonInterface wilson;
    wilson.build(WilsonBuildConfig(std::unordered_set<WGroup>{WGroup::B}, 81.0, 4.8, QCDOrder::NNLO));

    const std::vector<double> hadronic_scales = {2.0, 4.2, 4.8, 5.5};

    for (double mu_h : hadronic_scales) {
        // These setters update the corresponding scale parameters in the global parameter system.
        wilson.set_hadronic_scale(mu_h);

        const auto c9 = wilson.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL);
        const auto c10 = wilson.getFR(WGroup::B, WCoef::C10, QCDOrder::NNLO, ContributionType::TOTAL);

        std::cout << "mu_h = " << mu_h
                  << " -> C9(full) = " << c9
                  << ", C10(full) = " << c10 << "\n";
    }

    return 0;
}
