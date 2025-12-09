#ifndef INTERPRETED_PARAM_H
#define INTERPRETED_PARAM_H

#include "Include.h"

/**
 * @file InterpretedParam.h
 * @brief Lightweight key describing an interpreted LHA/SLHA parameter.
 *
 * InterpretedParam is used as a compact identifier for parameters that
 * have been read and interpreted from an LHA-like source. It combines:
 *   - a block name,
 *   - an LhaID code (possibly multi-index),
 *   - flags describing the model sector (SM/BSM) and whether the value
 *     is treated as complex.
 *
 * A std::hash specialization is provided so it can be used as a key in
 * unordered containers.
 */

/**
 * @struct InterpretedParam
 * @brief Identifies a single interpreted parameter from an LHA block.
 *
 * The struct does not store the actual numeric value, only the metadata
 * needed to uniquely identify which parameter is being referred to.
 */
struct InterpretedParam {
    /// Name of the LHA block (e.g. "SMINPUTS", "MASS", "B_Dlnu", ...).
    std::string block;

    /// LHA index (or multi-index) identifying the entry inside the block.
    LhaID code;

    /// True if this parameter belongs to the BSM part of the model.
    bool is_bsm;

    /// True if the parameter is interpreted as complex-valued.
    bool is_complex;

     /**
     * @brief Equality comparison operator.
     *
     * Two InterpretedParam objects are equal if all fields match exactly.
     */
    bool operator==(const InterpretedParam& other) const {
        return block == other.block &&
                code == other.code &&
                is_bsm == other.is_bsm &&
                is_complex == other.is_complex;
    }
};

/**
 * @brief Helper to combine hash values (boost-like hash_combine).
 *
 * This function is used by the std::hash specialization for
 * InterpretedParam to mix the component hashes into a single seed.
 *
 * @param seed  Hash seed to update.
 * @param value Hash value to combine into the seed.
 */
inline void hash_combine(std::size_t& seed, std::size_t value) noexcept {
    seed ^= value + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
}

/**
 * @brief std::hash specialization for InterpretedParam.
 *
 * The hash is computed by combining:
 *   - the hash of the block string,
 *   - the hash of the LhaID,
 *   - the hash of the is_bsm flag,
 *   - the hash of the is_complex flag.
 */
namespace std {
    template<>
    struct hash<InterpretedParam> {
        std::size_t operator()(const InterpretedParam& p) const noexcept {
            std::size_t seed = 0;
            hash_combine(seed, std::hash<std::string>{}(p.block));
            hash_combine(seed, hash<LhaID>{}(p.code));
            hash_combine(seed, std::hash<bool>{}(p.is_bsm));
            hash_combine(seed, std::hash<bool>{}(p.is_complex));
            return seed;
        }
    };
}


#endif