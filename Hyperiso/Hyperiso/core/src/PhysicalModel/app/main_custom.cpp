#include <iostream>
#include "WilsonBuilder.h"
#include "WilsonManager.h"
#include "GroupDefinition.h"
#include "CustomWilson.h"
#include "CustomWilsonGroup.h"

int main() {

    auto hyp = HyperisoMaster();
    Config config_hyp;

    config_hyp.model = Model::SM;
    config_hyp.mty_model_name = "ZPrime";
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h";
    hyp.init("lha/camilia.flha", config_hyp);

    GroupMapper::register_custom("NP_DEMO", {"np-demo"});
    const WGroupId GID = GroupMapper::id_of("NP_DEMO");

    WCoefMapper::register_custom("C_NP1");
    WCoefMapper::register_custom("C_NP2");
    const auto CNP1 = WCoefMapper::id_of("C_NP1");
    const auto CNP2 = WCoefMapper::id_of("C_NP2");
    WCoefMapper::set_external(CNP1, {9001, 1});
    WCoefMapper::set_external(CNP2, {9001, 2});

    GroupDefinition def;
    def.id = GID;

    def.sources[WilsonBasis::B_STANDARD][QCDOrder::LO].sources = {
        { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER } }
    };

    def.sources[WilsonBasis::B_STANDARD][QCDOrder::LO].func =
        [=](const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& match,
            const BlockSrc& /*src*/) -> std::unordered_map<WCoefId, scalar_t>
    {
        const auto& C_LO = match.at(QCDOrder::LO);

        const auto it1 = C_LO.find(CNP1);
        const auto it2 = C_LO.find(CNP2);
        scalar_t C1W = (it1==C_LO.end()) ? scalar_t{} : it1->second;
        scalar_t C2W = (it2==C_LO.end()) ? scalar_t{} : it2->second;

        // [ 0.85  0.10 ]
        // [ 0.05  0.90 ]
        scalar_t C1H = scalar_t(0.85)*C1W + scalar_t(0.10)*C2W;
        scalar_t C2H = scalar_t(0.05)*C1W + scalar_t(0.90)*C2W;

        return {
            { CNP1, C1H },
            { CNP2, C2H }
        };
    };

    def.setup[Model::SM].push_back([=](const BuildContext& ctx, CoefficientGroup& grp){
        const auto blk = GroupMapper::str(ctx.group_id, ScaleType::MATCHING);

        auto w1 = std::make_shared<CustomWilson>(LhaID(9001,1,0,(int)ctx.contrib), blk, QCDOrder::LO, ctx.contrib);
        w1->set_name("C_NP1");
        w1->set_contribution_type(ctx.contrib);
        w1->set_order_info(QCDOrder::LO, /*sources*/{},
            [](const ParamSrc&){ return scalar_t(1.234); },
            LhaID(9001,1,0,(int)ctx.contrib));

        auto w2 = std::make_shared<CustomWilson>(LhaID(9001,2,0,(int)ctx.contrib), blk, QCDOrder::LO, ctx.contrib);
        w2->set_name("C_NP2");
        w2->set_contribution_type(ctx.contrib);
        w2->set_order_info(QCDOrder::LO, /*sources*/{},
            [](const ParamSrc&){ return scalar_t(0.5); },
            LhaID(9001,2,0,(int)ctx.contrib));

        grp.insert({ w1->get_name(), w1 });
        grp.insert({ w2->get_name(), w2 });
    });

    GroupDefinitions::register_custom(def);

    WilsonBuildConfig cfg;
    cfg.groups = { GID };
    cfg.matching_scale = 80.0;  // mu_W
    cfg.hadronic_scale = 4.8;   // mu_h
    cfg.order = QCDOrder::LO;

    auto builder = std::make_shared<WilsonBuilder>(cfg);
    auto cm = builder->get_coefficient_manager();

    const std::string gname = GroupMapper::str(GID);
    auto CNP1_M = cm->getMatchingCoefficient(gname, "C_NP1", "LO", ContributionType::SM);
    auto CNP2_M = cm->getMatchingCoefficient(gname, "C_NP2", "LO", ContributionType::SM);
    auto CNP1_R = cm->getRunCoefficient(gname, "C_NP1", "LO", ContributionType::SM, WilsonBasis::B_STANDARD);
    auto CNP2_R = cm->getRunCoefficient(gname, "C_NP2", "LO", ContributionType::SM, WilsonBasis::B_STANDARD);

    std::cout << "C_NP1 LO @ mu_W = " << CNP1_M << "\n";
    std::cout << "C_NP2 LO @ mu_W = " << CNP2_M << "\n";
    std::cout << "C_NP1 LO @ mu_h = " << CNP1_R << "\n";
    std::cout << "C_NP2 LO @ mu_h = " << CNP2_R << "\n";

    return 0;
}
