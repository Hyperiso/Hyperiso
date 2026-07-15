#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>

#include "MemoryManager.h"
#include "ParameterProvider.h"
#include "Parameters.h"
#include "WilsonInterface.h"

namespace {

double real_value(const scalar_t& value) {
    return static_cast<complex_t>(value).real();
}

void print_complex(const std::string& label, const scalar_t& value) {
    const auto z = static_cast<complex_t>(value);
    std::cout << std::left << std::setw(28) << label << std::scientific
              << std::setprecision(16) << z.real()
              << (z.imag() >= 0.0 ? " +" : " ") << z.imag() << "i\n";
}

double optional_value(const ParameterProvider& p, const std::string& block, int id) {
    return p.exists(block, LhaID(id)) ? real_value(p(block, LhaID(id))) : 0.0;
}

} // namespace

int main(int argc, char** argv) {
    const std::string input = argc > 1 ? argv[1] : "lha/nmssm_smoke.slha";
    const double q_match = argc > 2 ? std::atof(argv[2]) : 81.0;
    const double q_low = argc > 3 ? std::atof(argv[3]) : q_match;

    try {
        auto hyp = HyperisoMaster();
        HyperisoConfig config_hyp;
        config_hyp.model = Model::SUSY;
        config_hyp.flags[ExternalFlag::IS_LHA_SPECTRUM] = true;
        hyp.init(input, config_hyp);

        const ParameterProvider bsm(ParameterType::BSM);
        const ParameterProvider pass(ParameterType::PASSTHROUGH);
        const auto nm = [&](int run_id, int ext_id) {
            if (bsm.exists("NMSSMRUN", LhaID(run_id))) return real_value(bsm("NMSSMRUN", LhaID(run_id)));
            if (pass.exists("EXTPAR", LhaID(ext_id))) return real_value(pass("EXTPAR", LhaID(ext_id)));
            return 0.0;
        };

        std::cout << std::scientific << std::setprecision(16)
                  << "HYPERISO_NMSSM_INPUT lambda=" << nm(1, 61)
                  << " kappa=" << nm(2, 62)
                  << " A_lambda=" << nm(3, 63)
                  << " mu_eff=" << nm(5, 65) << '\n'
                  << "HYPERISO_NMSSM_MASSES h3=" << optional_value(bsm, "MASS", 45)
                  << " a1=" << optional_value(bsm, "MASS", 36)
                  << " a2=" << optional_value(bsm, "MASS", 46) << '\n'
                  << "HYPERISO_YUKAWA_BLOCKS YU33="
                  << (bsm.exists("YU", LhaID(3, 3)) ? "present" : "derived")
                  << " YD33=" << (bsm.exists("YD", LhaID(3, 3)) ? "present" : "derived") << '\n'
                  << "HYPERISO_NEUTRALINO5 mass=" << optional_value(bsm, "MASS", 1000045)
                  << " matrix=" << (bsm.exists("NMNMIX", LhaID(1, 1)) ? "NMNMIX" : "NMIX")
                  << " row5="
                  << ((bsm.exists("NMNMIX", LhaID(5, 3)) && bsm.exists("NMNMIX", LhaID(5, 4))) ||
                      (bsm.exists("NMIX", LhaID(5, 3)) && bsm.exists("NMIX", LhaID(5, 4)))
                          ? "complete" : "incomplete") << '\n'
                  << "HYPERISO_SCALES q_match=" << q_match << " q_low=" << q_low << '\n';

        WilsonInterface wi;
        WilsonBuildConfig config(std::unordered_set<WGroup>{WGroup::BScalar, WGroup::B},
                                 q_match, q_low, QCDOrder::NLO);
        wi.build(config);

        const ParameterProvider wilson(ParameterType::WILSON);
        std::cout << "HYPERISO_EPSILON_SUSY eps0="
                  << real_value(wilson("EPSILON_SUSY", LhaID(0, 1)))
                  << " eps0p=" << real_value(wilson("EPSILON_SUSY", LhaID(0, 2)))
                  << " eps1p=" << real_value(wilson("EPSILON_SUSY", LhaID(1)))
                  << " eps2=" << real_value(wilson("EPSILON_SUSY", LhaID(2)))
                  << " epsb=" << real_value(wilson("EPSILON_SUSY", LhaID(3)))
                  << " epsbp=" << real_value(wilson("EPSILON_SUSY", LhaID(4))) << '\n';

        print_complex("HYPERISO_CQ1_LO_MATCH", wi.getM(WGroup::BScalar, WCoef::CQ1_MU, QCDOrder::LO, ContributionType::BSM));
        print_complex("HYPERISO_CQ1_NLO_MATCH", wi.getM(WGroup::BScalar, WCoef::CQ1_MU, QCDOrder::NLO, ContributionType::BSM));
        print_complex("HYPERISO_CQ2_LO_MATCH", wi.getM(WGroup::BScalar, WCoef::CQ2_MU, QCDOrder::LO, ContributionType::BSM));
        print_complex("HYPERISO_CQ2_NLO_MATCH", wi.getM(WGroup::BScalar, WCoef::CQ2_MU, QCDOrder::NLO, ContributionType::BSM));
        print_complex("HYPERISO_C7_LO_MATCH", wi.getM(WGroup::B, WCoef::C7, QCDOrder::LO, ContributionType::BSM));
        print_complex("HYPERISO_C7_NLO_MATCH", wi.getM(WGroup::B, WCoef::C7, QCDOrder::NLO, ContributionType::BSM));
        print_complex("HYPERISO_C8_LO_MATCH", wi.getM(WGroup::B, WCoef::C8, QCDOrder::LO, ContributionType::BSM));
        print_complex("HYPERISO_C8_NLO_MATCH", wi.getM(WGroup::B, WCoef::C8, QCDOrder::NLO, ContributionType::BSM));
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "HyperIso NMSSM test failed: " << e.what() << '\n';
        return 1;
    }
}
