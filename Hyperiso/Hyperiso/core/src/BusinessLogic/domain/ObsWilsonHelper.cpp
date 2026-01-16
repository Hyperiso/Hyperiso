#include "ObsWilsonHelper.h"

void ObsWilsonHelper::build(WilsonBuildConfig config, std::shared_ptr<IObsWilsonBuilder>& wil_builder, std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer) {
    config.groups = update_state(config.groups, iobs_wfreezer);
    if (config.groups.empty()) {
        return;
    }
    for (auto group : config.groups) {
        LOG_DEBUG("Building Wilson group", GroupMapper::str(group));
    }
    wil_builder->build(std::make_shared<WilsonBuildConfig>(config));
}

std::unordered_set<WGroupId> ObsWilsonHelper::get_all_groups(const std::unordered_set<WGroupId> &needed) {
    std::unordered_set<WGroupId> all_groups = get_keys(state);
    std::set_union(
        all_groups.begin(), all_groups.end(),
        needed.begin(), needed.end(),
        std::inserter(all_groups, all_groups.end())
    );
    return all_groups;
}

std::unordered_set<WGroupId> ObsWilsonHelper::update_state(const std::unordered_set<WGroupId> &needed, std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer) {
    std::unordered_set<WGroupId> to_build;
    for (auto group : get_all_groups(needed)) {
        if (!needed.contains(group)) {
            iobs_wfreezer->freeze(group);
            state[group] = true;
        } else if (!state.contains(group)) {
            to_build.emplace(group);
            state[group] = false;
        } else if (state[group]) {
            iobs_wfreezer->unfreeze(group);
            state[group] = false;
        }
    }

    return to_build;
}
