#include "CoefficientGroupBuilder.h"
#include "FileNameManager.h"

#include <filesystem>

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

    // Prepare one shared MARTY model for the active members whenever the
    // selected contribution is genuinely MARTY-backed.  For a pure SM MARTY
    // build, coefficients without a shipped MARTY template are omitted from
    // the analytical batch and are resolved by the registry's native-SM
    // fallback below.
    if (ctx.backend == Backend::Marty
        && ctx.adapters.marty_proxy) {
        if (ctx.contrib == ContributionType::BSM
            && ctx.model != Model::SM
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
                false, true, false
            );
        } else if (ctx.contrib == ContributionType::SM
                   && ctx.model == Model::SM) {
            std::vector<std::string> marty_members;
            marty_members.reserve(active_members.size());
            for (auto c : active_members) {
                // These SM branches are intentionally native even in a MARTY
                // session; do not spend analytical time batching a result the
                // registry will discard.
                const bool native_compat =
                    c == WCoef::C9 || c == WCoef::CP9 || c == WCoef::CP10
                    || def.id == GroupMapper::to_id(WGroup::BNuNu)
                    || def.id == GroupMapper::to_id(WGroup::KNuNu);
                if (native_compat) {
                    continue;
                }

                const std::string name = WCoefMapper::str(c);
                const auto files = FileNameManager::getInstance(name, "SM");
                const std::filesystem::path tmpl =
                    std::filesystem::path(files->getTemplateDir()) / (name + ".cpp");
                if (std::filesystem::is_regular_file(tmpl)) {
                    marty_members.push_back(name);
                }
            }
            if (!marty_members.empty()) {
                ctx.adapters.marty_proxy->prepare_group(
                    GroupMapper::str(def.id, ScaleType::MATCHING),
                    marty_members,
                    "SM", "SM", ctx.adapters.sm_path.string(),
                    false, false, false
                );
            }
        }
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

