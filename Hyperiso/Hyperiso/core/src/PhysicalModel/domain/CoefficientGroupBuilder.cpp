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

    std::vector<WCoef> active_members;
    active_members.reserve(def.members.size());
    for (auto c : def.members) {
        const auto id = WCoefMapper::to_id(c);
        if (ctx.requested_coefficients.empty()
            || ctx.requested_coefficients.contains(id)) {
            active_members.emplace_back(c);
        }
    }

    if (!ctx.requested_coefficients.empty() && active_members.empty()) {
        throw std::runtime_error(
            "Requested Wilson coefficient subset has no member in group '"
            + GroupMapper::str(def.id) + "'"
        );
    }

    std::vector<WCoefId> member_ids;
    member_ids.reserve(active_members.size());

    // For a genuine BSM MARTY contribution, prepare only the active members.
    // This lets matching-only diagnostics request C7/C8 without computing
    // unrelated one-loop four-fermion coefficients from the same group.
    if (ctx.backend == Backend::Marty
        && ctx.contrib == ContributionType::BSM
        && ctx.model != Model::SM
        && ctx.adapters.marty_proxy
        && ctx.adapters.marty_model_name
        && ctx.adapters.marty_model_path) {
        std::vector<std::string> marty_members;
        marty_members.reserve(active_members.size());
        for (auto c : active_members) marty_members.push_back(WCoefMapper::str(c));
        ctx.adapters.marty_proxy->prepare_group(
            GroupMapper::str(def.id, ScaleType::MATCHING),
            marty_members,
            ctx.adapters.marty_model_name->get(),
            ctx.adapters.marty_model_name->get(),
            ctx.adapters.marty_model_path->get().string(),
            false,
            true,
            false
        );
    }

    for (auto c : active_members) {
        auto coef = reg_.create(ctx, c);
        grp->insert({ WCoefMapper::str(c), std::move(coef) });
        member_ids.emplace_back(WCoefMapper::to_id(c));
    }

    grp->set_member_ids(std::move(member_ids));
    
    for (auto& hook : def.hooks_for(ctx.model)) {
        hook(ctx, *grp);
    }

    return grp;
}

