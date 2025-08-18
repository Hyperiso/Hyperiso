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
                {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)}, // x_t
                {ParameterType::WILSON, "EW_SCALE", 1}                  // Q_match
            },
            [](const auto& src)->scalar_t {
                double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)})->get_val();
                double Q  = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
                // démo: remplace par ta vraie formule
                return 0.123 + 0.5*std::log(Q) - 0.1*xt;
            },
            LhaID(999, 1001, 0, 0) // LHA ID LO
        );

        // NLO : sources + lambda + LHA ID
        Cfoo->set_order_info(
            QCDOrder::NLO,
            std::unordered_set<ParamId>{
                {ParameterType::WILSON, "WPARAM_MATCH_SM", 3},          // L
                {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)}, // x_t
                {ParameterType::WILSON, "EW_SCALE", 1}                  // Q_match
            },
            [](const auto& src)->scalar_t {
                double L  = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
                double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)})->get_val();
                return 0.456 + 0.25*L - 0.02*xt*xt;
            },
            LhaID(999, 1001, 1, 0) // LHA ID NLO
        );

        // 2) Créer un groupe custom B et y ajouter le coefficient
        auto Gcustom = std::make_shared<CustomCoefficientGroup>(WGroup::B, ContributionType::SM);
        Gcustom->add_coefficient(Cfoo);

        // 3) Déclarer les sources et le running pour la base standard à LO
        Gcustom->set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::LO,
            std::unordered_map<ParameterType, std::vector<std::string>>{
                {ParameterType::WILSON, {"U_MATRIX", "WPARAM_MATCH_SM", "WPARAM_RUN_SM"}},
                {ParameterType::SM,     {"MASS"}}
            },
            // Running “identité” (copie le matching de l’ordre dispo vers le hadronique)
            CustomCoefficientGroup::identity_running
        );

        // 4) Finaliser : compose_parameter sur chaque coef/ordre
        Gcustom->finalize(QCDOrder::NLO);

        // 5) Lire et afficher quelques valeurs
        // Matching mu_W
        auto cfoo_lo_match  = Gcustom->get_matching_coefficient("Cfoo", "LO",  ContributionType::TOTAL);
        auto cfoo_nlo_match = Gcustom->get_matching_coefficient("Cfoo", "NLO", ContributionType::TOTAL);

        // Running mu_h (base standard)
        auto cfoo_lo_run  = Gcustom->get_running_coefficient("Cfoo", "LO",  ContributionType::TOTAL, WilsonBasis::B_STANDARD);
        auto cfoo_nlo_run = Gcustom->get_running_coefficient("Cfoo", "NLO", ContributionType::TOTAL, WilsonBasis::B_STANDARD);

        std::cout << "=== Custom test ===\n";
        std::cout << "Cfoo (matching, LO)  = " << cfoo_lo_match  << "\n";
        std::cout << "Cfoo (matching, NLO) = " << cfoo_nlo_match << "\n";
        std::cout << "Cfoo (running,  LO)  = " << cfoo_lo_run    << "\n";
        std::cout << "Cfoo (running,  NLO) = " << cfoo_nlo_run   << "\n";

        // 6) Optionnel : afficher tout le groupe (utilise ton operator<<)
        std::cout << "\n--- Dump du groupe ---\n";
        std::cout << Gcustom << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
}
