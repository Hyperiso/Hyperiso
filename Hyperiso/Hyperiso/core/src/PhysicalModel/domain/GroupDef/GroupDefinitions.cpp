#include "GroupDefinition.h"

#include <algorithm>


std::vector<GroupDefinition::SetupHook> GroupDefinition::hooks_for(Model model) const {
    std::vector<SetupHook> hooks = common_setup;

    auto append_unique = [&hooks](const std::vector<SetupHook>& src) {
        for (const auto& hook : src) {
            const auto* target = hook.target<void(*)(const BuildContext&, CoefficientGroup&)>();
            const bool already_present = target && std::any_of(
                hooks.begin(), hooks.end(),
                [target](const SetupHook& existing) {
                    const auto* existing_target = existing.target<void(*)(const BuildContext&, CoefficientGroup&)>();
                    return existing_target && (*existing_target == *target);
                });
            if (!already_present) {
                hooks.push_back(hook);
            }
        }
    };

    if (model == Model::MARTY) {
        if (auto it = setup.find(Model::SM); it != setup.end()) {
            append_unique(it->second);
        }
    }

    if (auto it = setup.find(model); it != setup.end()) {
        append_unique(it->second);
    }

    return hooks;
}

namespace {
    static std::unordered_map<WGroupId, GroupDefinition> kCustomDefs;
}

namespace GroupDefinitions {
    void register_custom(const GroupDefinition& def) {
        kCustomDefs[def.id] = def;
    }
    bool has_custom(WGroupId id) {
        return kCustomDefs.find(id) != kCustomDefs.end();
    }
    const GroupDefinition& get_custom(WGroupId id) {
        auto it = kCustomDefs.find(id);
        if (it == kCustomDefs.end()) throw std::runtime_error("Unknown custom WGroupId");
        return it->second;
    }

    const GroupDefinition& get(WGroupId g) {
        if (auto en = GroupMapper::enum_of(g)) {
            switch (*en) {
                case WGroup::B:            return B();
                case WGroup::BPrime:       return BPrime();
                case WGroup::BScalar:      return BScalar();
                case WGroup::CC_bc:          return CC_bc();
                case WGroup::CC_bu:          return CC_bu();
                case WGroup::CC_cs:          return CC_cs();
                case WGroup::CC_cd:          return CC_cd();
                case WGroup::CC_su:          return CC_su();
                case WGroup::CC_du:          return CC_du();
                case WGroup::MESON_MIXING: return MesonMixing();
                case WGroup::K:            return K();
                default: break;
            }
            throw std::runtime_error("Unknown built-in WGroup");
        }

        return get_custom(g);
    }
}