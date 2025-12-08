#ifndef PARAM_ID_H
#define PARAM_ID_H

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
#include "LhaID.h"
#include "BlockName.h"

namespace fs = std::filesystem;

/**
 * @file ParamId.h
 * @brief Identifier for input parameters in LHA / SLHA-style structures.
 *
 * This header defines the ParamId struct, a lightweight key type used to
 * unambiguously identify a parameter in the context of LHA/SLHA blocks.
 *
 * A parameter is characterized by:
 *   - an optional high-level parameter type (ParameterType),
 *   - the name of the block it belongs to (BlockName),
 *   - an index or multi-index within that block (LhaID).
 *
 * ParamId is intended to be used as a key in maps and sets (including
 * unordered containers, thanks to the hash specialization).
 */

/**
 * @struct ParamId
 * @brief Composite identifier for a single parameter.
 *
 * ParamId ties together three pieces of information:
 *   - @ref type  : semantic category of the parameter (SM, FLAVOR, DECAY, ...),
 *   - @ref block : SLHA/LHA block name (with alias support),
 *   - @ref code  : index or multi-index inside the block (e.g. PDG code, (i,j), ...).
 *
 * This allows one to uniquely reference a given entry in a parameter file
 * or internal parameter store. The type field is optional; for some use
 * cases only the (block, code) pair is relevant.
 *
 * Comparison operators (== and <) are defined so that ParamId can be used
 * as a key in ordered containers. A hash specialization is also provided
 * to enable its use in unordered containers.
 */
struct ParamId {
    /**
     * @brief Optional high-level parameter category.
     *
     * This can be used to distinguish between, e.g. SM, flavor, decay or
     * Wilson parameter types (see ParameterType enum). If not set, only
     * the block and code fields are considered.
     */
    std::optional<ParameterType> type;

    /**
     * @brief Name of the block where the parameter is stored.
     *
     * This is typically an SLHA/LHA block name, possibly with aliases
     * (e.g. "SMINPUTS", "MASS", "DECAY", "FMASS", ...).
     */
    BlockName block;

    /**
     * @brief Index or multi-index of the parameter inside the block.
     *
     * This is represented by an LhaID, allowing both simple indices
     * (e.g. PDG codes) and multi-dimensional indices (i,j,...) to be
     * encoded in a uniform way.
     */
    LhaID code;

    /**
     * @brief Default constructor.
     *
     * Initializes a "null" parameter identifier with:
     *   - type unset,
     *   - block set to "NULL",
     *   - code set to 0.
     *
     * This can be used as a sentinel value.
     */
    ParamId();

    /**
     * @brief Constructs a ParamId from a block name and an index.
     *
     * The parameter type remains unset (std::nullopt). This is useful
     * when the (block, code) pair is sufficient to identify the parameter.
     *
     * @param block Block to which the parameter belongs.
     * @param code  Index or multi-index within the block.
     */
    ParamId(const BlockName& block, const LhaID& code);

    /**
     * @brief Constructs a ParamId from a type, block name and index.
     *
     * This constructor sets all three fields explicitly, providing the
     * most detailed identification of the parameter.
     *
     * @param type  High-level parameter category.
     * @param block Block to which the parameter belongs.
     * @param code  Index or multi-index within the block.
     */
    ParamId(ParameterType type, const BlockName& block, const LhaID& code);

    /**
     * @brief Sets or overwrites the parameter type.
     *
     * This can be useful when constructing a ParamId in several steps
     * or when inferring the parameter type after the fact.
     *
     * @param type New parameter type to assign.
     */
    void set_parameter_type(ParameterType type);

    /**
     * @brief Equality comparison between two ParamId objects.
     *
     * Two ParamId objects are equal if and only if all three components
     * (type, block, code) are equal. For the type field, equality of
     * std::optional<ParameterType> is used, so both must either be unset
     * or hold the same value.
     *
     * @param lhs Left-hand side.
     * @param rhs Right-hand side.
     * @return true if the identifiers are identical, false otherwise.
     */
    bool friend operator==(const ParamId& lhs, const ParamId& rhs);

    /**
     * @brief Strict ordering operator for ParamId.
     *
     * The ordering is defined lexicographically on (type, block, code):
     *   - First compare @ref type; if they differ, the result is based
     *     on the ordering of std::optional<ParameterType>.
     *   - If @ref type are equal, compare @ref block using BlockName::operator<.
     *   - If both type and block are equal, compare @ref code using LhaID::operator<.
     *
     * This allows ParamId to be used as a key in ordered containers such
     * as std::map or std::set.
     *
     * @param other Right-hand side.
     * @return true if @c *this is lexicographically smaller than @p other.
     */

    bool operator<(const ParamId& other) const;
};

/**
 * @brief Hash specialization for ParamId.
 *
 * This allows ParamId to be used as a key in unordered associative
 * containers such as std::unordered_map or std::unordered_set.
 *
 * The hash is computed by combining the hashes of:
 *   - the BlockName (block),
 *   - the LhaID (code),
 *   - the optional ParameterType (type, if present).
 */
namespace std {
    template <>
    struct hash<ParamId> {
        /**
         * @brief Computes a hash value for a ParamId.
         *
         * The three components are hashed and combined via simple XOR and
         * bit-shifting. This is sufficient for typical usage but does not
         * claim cryptographic properties.
         *
         * @param p ParamId to hash.
         * @return Combined hash value.
         */
        std::size_t operator()(const ParamId& p) const noexcept {
            std::size_t h1 = std::hash<BlockName>{}(p.block);
            std::size_t h2 = std::hash<LhaID>{}(p.code);
            std::size_t h3 = p.type ? std::hash<int>{}(static_cast<int>(*p.type)) : 0;
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

/**
 * @brief Stream output operator for ParamId.
 *
 * Prints the parameter identifier in the form:
 *   `<block>:<code>`
 *
 * where both parts use their respective stream operators. The type
 * field is not printed.
 *
 * Example:
 * @code
 * ParamId p(ParameterType::SM, "MASS", LhaID(5));
 * std::cout << p; // might print: MASS:5
 * @endcode
 *
 * @param os  Output stream.
 * @param pid ParamId to print.
 * @return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const ParamId& pid) {
    os << pid.block << ":" << pid.code;
    return os;
};

#endif