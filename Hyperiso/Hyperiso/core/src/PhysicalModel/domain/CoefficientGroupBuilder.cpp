#include "CoefficientGroupBuilder.h"

/**
 * @file CoefficientGroupBuilder.cpp
 * @brief Implementation of @ref CoefficientGroupBuilder.
 *
 * The implementation follows the build steps documented in the header.
 * In particular, it performs placeholder substitution for sources and
 * executes model-specific setup hooks after all members are created.
 */

std::shared_ptr<CoefficientGroup> CoefficientGroupBuilder::build(const BuildContext& ctx) const {
    auto def = GroupDefinitions::get(ctx.group_id);

    auto grp = std::make_shared<GenericCoefficientGroup>(ctx.adapters);
    grp->set_group_id(def.id);
    grp->set_wilson_type(ctx.contrib);

    std::string matching_block =
        (ctx.group_name.empty())
            ? GroupMapper::str(def.id, ScaleType::MATCHING)
            : ctx.group_name;
         
    grp->set_matching_storage_block(matching_block);

    for (const auto& [basis, per_order] : def.sources) {
        std::map<QCDOrder, CoefficientGroupSources> m;
        for (const auto& [ord, s] : per_order) {
            auto s2 = s;
            for (auto& [ptype, names] : s2.sources)
                for (auto& n : names)
                    if (n == MATCHING_BLOCK_PLACEHOLDER) n = matching_block;
            m.emplace(ord, std::move(s2));
        }
        grp->add_sources(basis, m);
    }

    for (auto c : def.members) {
        auto coef = reg_.create(ctx, c);
        grp->insert({ WCoefMapper::str(c), std::move(coef) });
    }
    
    if (auto it = def.setup.find(ctx.model); it != def.setup.end()) {
        for (auto& hook : it->second) hook(ctx, *grp);
    }

    return grp;
}

