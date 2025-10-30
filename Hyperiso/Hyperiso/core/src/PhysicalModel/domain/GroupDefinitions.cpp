#include "GroupDefinition.h"
#include <unordered_map>
#include <stdexcept>

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
                case WGroup::BCC_bc:          return BCC_bc();
                case WGroup::BCC_bu:          return BCC_bu();
                case WGroup::BCC_cs:          return BCC_cs();
                case WGroup::BCC_cd:          return BCC_cd();
                case WGroup::BCC_su:          return BCC_su();
                case WGroup::MESON_MIXING: return MesonMixing();
                case WGroup::K:            return K();
                default: break;
            }
            throw std::runtime_error("Unknown built-in WGroup");
        }

        return get_custom(g);
    }
}