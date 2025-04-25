#include "ObsWilsonHelper.h"

void ObsWilsonHelper::build(WilsonConfig config) {
    config.groups = update_state(config.groups);
    WilsonAdapter().build(config);
}

std::unordered_set<WGroup> ObsWilsonHelper::get_all_groups(const std::unordered_set<WGroup> &needed) {
    std::unordered_set<WGroup> all_groups = get_keys(state);
    std::set_union(
        all_groups.begin(), all_groups.end(),
        needed.begin(), needed.end(),
        std::inserter(all_groups, all_groups.end())
    );
    return all_groups;
}

std::unordered_set<WGroup> ObsWilsonHelper::update_state(const std::unordered_set<WGroup> &needed) {
    std::unordered_set<WGroup> to_build;
    for (auto group : get_all_groups(needed)) {
        if (!needed.contains(group)) {
            WilsonFreezer().freeze(group);
            state[group] = true;
        } else if (!state.contains(group)) {
            to_build.emplace(group);
            state[group] = false;
        } else if (state[group]) {
            WilsonFreezer().unfreeze(group);
            state[group] = false;
        }
    }

    return to_build;
}
