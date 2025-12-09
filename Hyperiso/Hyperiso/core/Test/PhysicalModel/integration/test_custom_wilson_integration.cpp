#include <cassert>
#include <iostream>
#include <unordered_map>

#include "CustomWilson.h"
#include "CustomWilsonGroup.h"
#include "GroupDefinition.h"
#include "GroupMapper.h"
#include "wcoef_ids.hpp"
#include "SourcesView.h"
#include "Block.h"

#include "stubs.hpp"

int main() {
    std::cout << "== CustomWilson / CustomCoefficientGroup INTEGRATION ==\n";

    auto ctx = make_ctx(
        Model::SUSY,
        Backend::Builtin,
        ContributionType::BSM,
        GroupMapper::to_id(WGroup::B)
    );

    WGroupId custom_id = ctx.group_id; 
    CustomCoefficientGroup grp(ctx.adapters, custom_id, "MY_CUSTOM_BSM_GROUP", ContributionType::BSM);

    CustomWilson c_lo(LhaID(1001), "W_CUSTOM");
    c_lo.with_type(ContributionType::BSM)
        .with_storage_block("W_CUSTOM_BLOCK");

    CustomWilson c_nlo(LhaID(1002), "W_CUSTOM");
    c_nlo.with_type(ContributionType::BSM)
         .with_storage_block("W_CUSTOM_BLOCK");

    ParamId pid{ParameterType::SM, "MASS", LhaID(5)};
    LhaID lhaid_NLO(777);
    c_nlo.set_order_info(
        QCDOrder::NLO,
        {pid},
        [](const ParamSrc&) { return 0.0; },
        lhaid_NLO
    );

    grp.add_coefficient(std::make_shared<CustomWilson>(c_lo));
    grp.add_coefficient(std::make_shared<CustomWilson>(c_nlo));

    assert(grp.size() == 2);

    std::unordered_map<ParameterType, std::vector<std::string>> source_names{
        {ParameterType::SM,     {"MASS", "GAUGE"}},
        {ParameterType::WILSON, {"W_CUSTOM_BLOCK"}}
    };

    grp.set_basis_order_sources_and_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::LO,
        source_names,
        &CustomCoefficientGroup::identity_running
    );

    grp.set_basis_order_sources_and_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::NLO,
        source_names,
        &CustomCoefficientGroup::identity_running
    );

    grp.finalize(QCDOrder::NLO);
    assert(grp.get_order() == QCDOrder::NLO);

    auto bases = grp.get_bases();
    assert(bases.find(WilsonBasis::B_STANDARD) != bases.end());
    auto srcs_LO = grp.get_sources(QCDOrder::LO, WilsonBasis::B_STANDARD);
    assert(!srcs_LO.empty());

    {
        std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> coef_matching;
        auto idC1 = WCoefMapper::to_id(WCoef::C1);

        coef_matching[QCDOrder::LO][idC1]  = 1.0;
        coef_matching[QCDOrder::NLO][idC1] = 5.0;

        std::unordered_map<std::string, std::shared_ptr<Block>> dummy_blocks;
        BlockSrc src(dummy_blocks, "DUMMY");

        auto run = CustomCoefficientGroup::identity_running(coef_matching, src);
        assert(run.size() == 1);
        assert(run.at(idC1) == scalar_t(5.0));
    }

    std::cout << " CustomWilson / CustomCoefficientGroup integration suite passed.\n";
    return 0;
}
