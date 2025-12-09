#include <cassert>
#include <iostream>
#include <unordered_set>
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
    std::cout << "== CustomWilson / CustomCoefficientGroup UNIT ==\n";

    {
        LhaID base_name(10); 
        std::string storage = "W_CUSTOM";

        CustomWilson cw(base_name, storage);

        cw.with_type(ContributionType::BSM)
          .with_storage_block("W_NEW_BLOCK");

        assert(cw.get_type() == ContributionType::BSM);
        assert(cw.get_storage_block() == "W_NEW_BLOCK");

        ParamId pid{ParameterType::SM, "MASS", LhaID(5)};
        std::unordered_set<ParamId> sources = {pid};

        LhaID lhaid_NLO(42);

        cw.set_order_info(
            QCDOrder::NLO,
            std::move(sources),
            [](const ParamSrc&) { return 1.234; },
            lhaid_NLO
        );

        auto srcs_NLO = cw.get_sources(QCDOrder::NLO);
        assert(srcs_NLO.size() == 1);
        assert(srcs_NLO.find(pid) != srcs_NLO.end());

        auto lhaNLO = cw.get_lhaid(QCDOrder::NLO);
        assert(lhaNLO.to_string() == lhaid_NLO.to_string());

        auto cloned = cw.clone();
        assert(cloned);
        assert(cloned->get_type() == ContributionType::BSM);
        assert(cloned->get_storage_block() == "W_NEW_BLOCK");

        auto srcs_NLO_clone = cloned->get_sources(QCDOrder::NLO);
        assert(srcs_NLO_clone.size() == 1);
        assert(srcs_NLO_clone.find(pid) != srcs_NLO_clone.end());

        auto lhaNLO_clone = cloned->get_lhaid(QCDOrder::NLO);
        assert(lhaNLO_clone.to_string() == lhaid_NLO.to_string());
    }

    {
        auto ctx = make_ctx(
            Model::SM,
            Backend::Builtin,
            ContributionType::SM,
            GroupMapper::to_id(WGroup::B) 
        );

        WGroupId gid_custom{};
        CustomCoefficientGroup grp(ctx.adapters, gid_custom, "MY_CUSTOM_GROUP", ContributionType::SM);

        CustomWilson cw(LhaID(100), "W_CUSTOM");
        cw.with_type(ContributionType::SM)
          .with_storage_block("W_CUSTOM_BLOCK");
        auto cw_ptr = std::make_shared<CustomWilson>(cw);

        grp.add_coefficient(cw_ptr);
        assert(grp.size() == 1);
        assert(grp.find(cw_ptr->get_name()) != grp.end());

        std::unordered_map<ParameterType, std::vector<std::string>> source_names{
            {ParameterType::SM,     {"MASS", "GAUGE"}},
            {ParameterType::WILSON, {"W_BLOCK"}}
        };

        grp.set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::LO,
            source_names,
            &CustomCoefficientGroup::identity_running
        );

        auto bases = grp.get_bases();
        assert(bases.find(WilsonBasis::B_STANDARD) != bases.end());

        auto srcs_LO = grp.get_sources(QCDOrder::LO, WilsonBasis::B_STANDARD);
        assert(srcs_LO.size() == 2);
        assert(srcs_LO.find(ParameterType::SM)     != srcs_LO.end());
        assert(srcs_LO.find(ParameterType::WILSON) != srcs_LO.end());

        grp.set_basis_order_sources_and_running(
            WilsonBasis::B_STANDARD,
            QCDOrder::NLO,
            source_names,
            &CustomCoefficientGroup::identity_running
        );

        assert(grp.get_order() == QCDOrder::NLO);

        std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>> coef_matching;
        auto idC1 = WCoefMapper::to_id(WCoef::C1);
        coef_matching[QCDOrder::LO][idC1]  = 1.0;
        coef_matching[QCDOrder::NLO][idC1] = 2.0;

        std::unordered_map<std::string, std::shared_ptr<Block>> dummy_blocks;
        BlockSrc src(dummy_blocks, "DUMMY");

        auto out = CustomCoefficientGroup::identity_running(coef_matching, src);
        assert(out.size() == 1);
        assert(out.at(idC1) == scalar_t(2.0));

        coef_matching[QCDOrder::NNLO][idC1] = 3.0;
        auto out2 = CustomCoefficientGroup::identity_running(coef_matching, src);
        assert(out2.at(idC1) == scalar_t(3.0));
    }

    std::cout << " CustomWilson / CustomCoefficientGroup unit suite passed.\n";
    return 0;
}
