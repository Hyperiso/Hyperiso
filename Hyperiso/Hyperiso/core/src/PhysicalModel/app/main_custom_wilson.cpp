// main.cpp
#include <iostream>
#include <memory>

#include "CustomWilson.h"
#include "CustomWilsonGroup.h"

#include "GroupMapper.h"
#include "QCDOrderMapper.h"
#include "WCoeffMapper.h"

// Optionnel: si tu as des logs
// #include "Utils.h"

int main() {

    auto hyp = HyperisoMaster();
    Config config_hyp;
    config_hyp.flags[ExternalFlag::USE_MARTY] = false;
    config_hyp.model = Model::SM;
    config_hyp.mty_model_name = "ZPrime";
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h";

    hyp.init("lha/camilia.flha", config_hyp); // Initialize program manager with LHA file


    try {
        LhaID a = {42,42};
        
        std::cout << "=== Custom Wilson Coefficient Test ===" << std::endl;
        auto Cfoo = std::make_shared<CustomWilson>(
            a,
            GroupMapper::str(WGroup::B) + std::string("_MATCH"),
            QCDOrder::NLO,
            ContributionType::SM
        );
        std::cout << "=== Creation done ===" << std::endl;

        Cfoo->set_order_info(
            QCDOrder::LO,
            std::unordered_set<ParamId>{
            },
            [](const auto& src)->scalar_t {
                double xt = 1;
                double Q = 4.18;
                return 0.123 + 0.5*std::log(Q) - 0.1*xt;
            },
            LhaID(999, 1001, 0, 0) 
        );

        Cfoo->set_order_info(
            QCDOrder::NLO,
            std::unordered_set<ParamId>{
            },
            [](const auto& src)->scalar_t {

                double L = 2;
                double xt = 1;
                return 0.456 + 0.25*L - 0.02*xt*xt;
            },
            LhaID(999, 1001, 1, 0)
        );

        auto Gcustom = std::make_shared<CustomCoefficientGroup>(WGroup::B, ContributionType::SM);
        Gcustom->add_coefficient(Cfoo);

        Gcustom->set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::LO,
            std::unordered_map<ParameterType, std::vector<std::string>>{
                {ParameterType::WILSON, {"U_MATRIX", "WPARAM_MATCH_SM", "WPARAM_RUN_SM"}},
                {ParameterType::SM,     {"MASS"}}
            },
            CustomCoefficientGroup::identity_running
        );

        Gcustom->finalize(QCDOrder::NLO);


        auto cfoo_lo_match  = Gcustom->get_matching_coefficient("Cfoo", "LO",  ContributionType::TOTAL);
        auto cfoo_nlo_match = Gcustom->get_matching_coefficient("Cfoo", "NLO", ContributionType::TOTAL);

        auto cfoo_lo_run  = Gcustom->get_running_coefficient("Cfoo", "LO",  ContributionType::TOTAL, WilsonBasis::B_STANDARD);
        auto cfoo_nlo_run = Gcustom->get_running_coefficient("Cfoo", "NLO", ContributionType::TOTAL, WilsonBasis::B_STANDARD);

        std::cout << "=== Custom test ===\n";
        std::cout << "Cfoo (matching, LO)  = " << cfoo_lo_match  << "\n";
        std::cout << "Cfoo (matching, NLO) = " << cfoo_nlo_match << "\n";
        std::cout << "Cfoo (running,  LO)  = " << cfoo_lo_run    << "\n";
        std::cout << "Cfoo (running,  NLO) = " << cfoo_nlo_run   << "\n";

        std::cout << "\n--- Dump du groupe ---\n";
        std::cout << Gcustom << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
}
