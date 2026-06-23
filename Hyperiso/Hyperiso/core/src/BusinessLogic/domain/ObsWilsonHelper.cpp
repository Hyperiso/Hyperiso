#include "ObsWilsonHelper.h"

void ObsWilsonHelper::clear() {
    state.clear();
}

void ObsWilsonHelper::build(WilsonBuildConfig config,
                            std::shared_ptr<IObsWilsonBuilder>& wil_builder,
                            std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer) {
    config.groups = update_state(config, std::move(iobs_wfreezer));
    if (config.groups.empty()) {
        return;
    }
    for (auto group : config.groups) {
        LOG_DEBUG("Building Wilson group", GroupMapper::str(group));
    }
    wil_builder->build(std::make_shared<WilsonBuildConfig>(config));
}

void ObsWilsonHelper::mark_built(WGroupId group, const WilsonBuildConfig& config, bool frozen) {
    state[group] = GroupState{frozen, make_signature(config)};
}

void ObsWilsonHelper::mark_built(WGroupId group,
                                 double matching_scale,
                                 double hadronic_scale,
                                 QCDOrder order,
                                 bool frozen) {
    WilsonBuildConfig config;
    config.groups = {group};
    config.matching_scale = matching_scale;
    config.hadronic_scale = hadronic_scale;
    config.order = order;
    mark_built(group, config, frozen);
}

std::unordered_set<WGroupId> ObsWilsonHelper::get_all_groups(const std::unordered_set<WGroupId>& needed) const {
    std::unordered_set<WGroupId> all_groups;
    all_groups.reserve(state.size() + needed.size());

    for (const auto& [group, _] : state) {
        all_groups.emplace(group);
    }
    all_groups.insert(needed.begin(), needed.end());

    return all_groups;
}

std::unordered_set<WGroupId> ObsWilsonHelper::update_state(const WilsonBuildConfig& config,
                                                           std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer) {
    std::unordered_set<WGroupId> to_build;
    const auto signature = make_signature(config);

    for (auto group : get_all_groups(config.groups)) {
        if (!config.groups.contains(group)) {
            auto it = state.find(group);
            if (it != state.end() && !it->second.frozen && iobs_wfreezer) {
                iobs_wfreezer->freeze(group);
                it->second.frozen = true;
            }
            continue;
        }

        auto it = state.find(group);
        if (it == state.end()) {
            to_build.emplace(group);
            state[group] = GroupState{false, signature};
            continue;
        }

        if (!it->second.signature.matches(config)) {
            to_build.emplace(group);
            it->second = GroupState{false, signature};
            continue;
        }

        if (it->second.frozen) {
            if (iobs_wfreezer) {
                iobs_wfreezer->unfreeze(group);
            }
            it->second.frozen = false;
        }
    }

    return to_build;
}

ObsWilsonHelper::BuildSignature ObsWilsonHelper::make_signature(const WilsonBuildConfig& config) {
    return BuildSignature{config.matching_scale, config.hadronic_scale, config.order};
}
