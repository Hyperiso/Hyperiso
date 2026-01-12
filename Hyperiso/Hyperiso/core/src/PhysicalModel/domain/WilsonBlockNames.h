#ifndef WILSON_BLOCK_NAME
#define WILSON_BLOCK_NAME

#include <string>

#include "wgroup_ids.hpp"
#include "scaletype_ids.hpp"

/**
 * @file WilsonBlockNames.h
 * @brief Centralizes canonical SLHA/LHA block naming conventions for Wilson blocks.
 *
 * This header defines small helpers that standardize how Wilson coefficient blocks
 * are named in the framework.
 *
 * The naming is based on GroupMapper + ScaleType conventions:
 *  - Matching blocks:   GroupMapper::str(gid, ScaleType::MATCHING)
 *  - SM-only matching:  "<matching>_SM"
 *  - BSM-only matching: "<matching>_BSM"
 *
 * It also exposes the fixed name for the FWCOEF block.
 *
 * @see GroupMapper
 * @see ScaleType
 * @see WGroupId
 */
struct WilsonBlockNames {
    /**
     * @brief Returns the canonical matching block name for a Wilson group.
     * @param gid Wilson group identifier.
     * @return Matching block name.
     */
    static std::string matching(WGroupId gid) {
        return GroupMapper::str(gid, ScaleType::MATCHING);
    }

    /**
     * @brief Returns the SM-only matching block name for a Wilson group.
     * @param gid Wilson group identifier.
     * @return "<matching>_SM".
     */
    static std::string sm_matching(WGroupId gid) {
        return matching(gid) + "_SM";
    }

    /**
     * @brief Returns the BSM-only matching block name for a Wilson group.
     * @param gid Wilson group identifier.
     * @return "<matching>_BSM".
     */
    static std::string bsm_matching(WGroupId gid) {
        return matching(gid) + "_BSM";
    }

    /**
     * @brief Returns the fixed FWCOEF block name.
     * @return "FWCOEF".
     */
    static constexpr const char* fwcoef() { return "FWCOEF"; }
};

#endif