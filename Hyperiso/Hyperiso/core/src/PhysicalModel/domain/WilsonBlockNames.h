#pragma once
#include "wgroup_ids.hpp"
#include "scaletype_ids.hpp"

struct WilsonBlockNames {
    static std::string matching(WGroupId gid) {
        return GroupMapper::str(gid, ScaleType::MATCHING);
    }

    static std::string sm_matching(WGroupId gid) {
        return matching(gid) + "_SM";
    }

    static std::string bsm_matching(WGroupId gid) {
        return matching(gid) + "_BSM";
    }

    static constexpr const char* fwcoef() { return "FWCOEF"; }
};
