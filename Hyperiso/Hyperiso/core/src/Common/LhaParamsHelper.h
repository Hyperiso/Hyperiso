#ifndef LHA_PARAMS_HELPER_H
#define LHA_PARAMS_HELPER_H

#include <string>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <initializer_list>
#include <map>
#include <vector>
#include <set>
#include <ranges>
#include <concepts>
#include <variant>
#include "Logger.h"
#include "Utils.h"
#include "GeneralEnum.h"
#include "BlockName.h"

namespace fs = std::filesystem;

/**
 * @file LhaParamsHelper.h
 * @brief Helper utilities for minimal LHAPDF / SLHA block content.
 *
 * This header defines the LhaParamsHelper class, which provides access
 * to the minimal set of entries that are required (or expected) to be
 * present in certain LHA/SLHA-like blocks.
 *
 * The minimal content is stored as:
 *   - a static map from BlockName to a vector of index-vectors,
 *   - where each entry corresponds to an identifier or multi-index
 *     (encoded as std::vector<long>) that should appear in the block.
 *
 * For example, for hadronic blocks like "FMASS" or "FLIFE" this may
 * encode which PDG codes must at least be provided.
 */

/**
 * @class LhaParamsHelper
 * @brief Helper for retrieving minimal required content of LHA/SLHA blocks.
 *
 * LhaParamsHelper exposes a single public static function:
 *   - get_minimal_content(const BlockName&),
 * which looks up the block in an internal table and returns a list of
 * required indices.
 *
 * The associated static map minimal_blocks is defined in the
 * implementation file and currently contains entries for a subset of
 * blocks such as:
 *   - "FMASS"  : minimal set of flavor masses,
 *   - "FLIFE"  : minimal set of lifetimes,
 *   - "FCONST" : minimal set of decay constants.
 *
 * This helper is typically used when validating an input parameter
 * file or when constructing a minimal template to be filled by the user.
 */
class LhaParamsHelper {
public:
    /**
     * @brief Returns the minimal required content for a given block.
     *
     * The supplied @p block_name is matched against the keys of the
     * internal minimal_blocks map (using BlockName comparison semantics,
     * including alias support). If an entry is found, the corresponding
     * list of index-vectors is returned.
     *
     * Each inner std::vector<long> represents an identifier or index
     * within the block. For example:
     *   - For "FMASS" and "FLIFE" the vectors may contain PDG codes,
     *     e.g. {211}, {321}, {511}, etc.
     *   - For "FCONST" the vectors may encode multi-indices, e.g.
     *     {511, 1}, {323, 2}, etc.
     *
     * If the block is not present in minimal_blocks, an error is
     * reported via LOG_ERROR in the implementation.
     *
     * @param block_name Block for which the minimal content is requested.
     * @return A vector of index-vectors describing the minimal content.
     */
    static std::vector<std::vector<long>> get_minimal_content(const BlockName& block_name);

private:
    /**
     * @brief Static lookup table of minimal required entries per block.
     *
     * The map keys are BlockName objects (with potential aliases),
     * and each mapped value is a list of minimal indices that should
     * appear in the corresponding block.
     *
     * This table is defined and initialized in LhaParamsHelper.cpp.
     * Example (non-exhaustive):
     * @code
     * {
     *   "FMASS"  : {{211}, {321}, {323}, {411}, ...},
     *   "FLIFE"  : {{211}, {321}, {323}, {411}, ...},
     *   "FCONST" : {{511,1}, {521,1}, {531,1}, {323,1}, {323,2}}
     * }
     * @endcode
     */
    static const std::map<BlockName, std::vector<std::vector<long>>> minimal_blocks;
};

#endif 