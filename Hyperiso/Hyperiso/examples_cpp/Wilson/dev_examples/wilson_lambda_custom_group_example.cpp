#include <iostream>
#include <unordered_set>

#include "HyperisoMaster.h"
#include "WilsonInterface.h"
#include "mapper_hub.hpp"

int main() {

    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);

    // Exemple dev : créer un groupe Wilson et deux coefficients entièrement dynamiques.
    // Les noms sont enregistrés à l'exécution : ils ne doivent pas exister dans les enums WGroup/WCoef.
    GroupMapper::register_custom("DEV_LAMBDA_GROUP", {"dev-lambda"});
    WCoefMapper::register_custom("C_DEV_1", {"cdev1"}, std::pair<int,int>{990001, 1});
    WCoefMapper::register_custom("C_DEV_2", {"cdev2"}, std::pair<int,int>{990001, 2});

    WGroupId group = GroupMapper::id_of("DEV_LAMBDA_GROUP");
    WCoefId c1 = WCoefMapper::id_of("C_DEV_1");
    WCoefId c2 = WCoefMapper::id_of("C_DEV_2");

    // Chaque coefficient reçoit ses lambdas de matching.
    // Ici on ne déclare aucune dépendance de paramètre, donc la lambda ignore ParamSrc.
    CustomWilsonCoefficientConfig coef1(c1);
    coef1.set_matching(
        QCDOrder::LO,
        {},
        [](const ParamSrc&) { return 1.25; },
        ContributionType::SM
    );

    CustomWilsonCoefficientConfig coef2(c2);
    coef2.set_matching(
        QCDOrder::LO,
        {},
        [](const ParamSrc&) { return -0.40; },
        ContributionType::SM
    );

    CustomWilsonGroupConfig cfg(group);
    cfg.matching_scale = 160.0;
    cfg.hadronic_scale = 4.8;
    cfg.order = QCDOrder::LO;
    cfg.contribution = ContributionType::SM;
    cfg.add_coefficient(coef1).add_coefficient(coef2);

    // Optionnel : définir l'évolution hadronique du groupe avec une lambda.
    // Si tu ne fournis rien, l'API installe une évolution identité en base B_STANDARD.
    cfg.set_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::LO,
        {},
        [c1, c2](const auto& matching, const BlockSrc&) {
            std::unordered_map<WCoefId, scalar_t> out;
            const auto& lo = matching.at(QCDOrder::LO);
            out[c1] = lo.at(c1);        // C_DEV_1(mu_b) = C_DEV_1(mu_W)
            out[c2] = 2.0 * lo.at(c2);  // exemple arbitraire de running custom
            return out;
        }
    );

    WilsonInterface W;
    W.addCustomWilsonGroup(cfg);

    std::cout << "C_DEV_1 matching LO = "
              << W.getM(group, c1, QCDOrder::LO, ContributionType::SM) << "\n";
    std::cout << "C_DEV_2 running LO = "
              << W.getR(group, c2, QCDOrder::LO, ContributionType::SM) << "\n";

    return 0;
}
