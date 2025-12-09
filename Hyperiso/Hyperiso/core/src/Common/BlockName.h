#ifndef BLOCK_NAME_H
#define BLOCK_NAME_H

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

namespace fs = std::filesystem;

/**
 * @file BlockName.h
 * @brief Representation of SLHA / LHA block names with support for aliases.
 *
 * This header defines the BlockName class, which encapsulates one or more
 * textual names referring to the same logical block (e.g. SLHA BLOCK names).
 *
 * The main features are:
 *   - storage of a set of aliases for a given block,
 *   - equality based on overlapping aliases (two BlockName are equal if they
 *     share at least one alias),
 *   - convenient construction from strings, initializer lists or sets,
 *   - implicit conversion to std::string for trivial (single-alias) blocks,
 *   - case normalization (conversion to upper case),
 *   - hash support for use in unordered containers.
 */

/**
 * @class BlockName
 * @brief Block identifier with alias support.
 *
 * A BlockName internally stores a set of strings, each string being an
 * alias that can refer to the same physical block. This is useful when
 * supporting multiple conventions or legacy naming schemes (e.g. "SMINPUTS"
 * vs "SMINPUT" or different capitalization).
 *
 * Key design points:
 *   - The internal representation is an unordered_set of strings.
 *   - Two BlockName instances are considered equal if they share at least
 *     one alias (see operator==(const BlockName&)).
 *   - Implicit conversion to std::string returns an arbitrary alias
 *     (the first in the internal set) and issues a warning if more than
 *     one alias is present.
 *   - The to_upper() method can be used to normalize all aliases to
 *     upper case, which is often convenient for SLHA/LHA parsing.
 */
class BlockName {
private:
    /**
     * @brief Set of all aliases for the block.
     *
     * Each entry is a string representation of a name that can identify
     * this block. The container is unordered; no ordering guarantees are given.
     */
    std::unordered_set<std::string> block_names;

public:
    /**
     * @brief Default constructor.
     *
     * Creates an empty BlockName without any aliases.
     */
    BlockName() = default;

    /**
     * @brief Constructs a BlockName from a single std::string.
     *
     * The provided string is inserted as the only alias of the block.
     *
     * @param name Initial block name / alias.
     */
    BlockName(const std::string& name);

    /**
     * @brief Constructs a BlockName from a C-style string.
     *
     * The provided C-string is converted to std::string and stored
     * as the only alias of the block.
     *
     * @param name Initial block name / alias (C-string).
     */
    BlockName(const char* name);
    
    /**
     * @brief Constructs a BlockName from an initializer list of aliases.
     *
     * All strings in @p names are inserted into the internal alias set.
     *
     * Example:
     * @code
     * BlockName b{"SMINPUTS", "sminputs"};
     * @endcode
     *
     * @param names List of aliases to associate with the block.
     */
    BlockName(std::initializer_list<std::string> names);

    /**
     * @brief Constructs a BlockName from an existing set of aliases.
     *
     * The provided set is copied into the internal representation.
     *
     * @param names Set of aliases for the block.
     */
    BlockName(const std::unordered_set<std::string>& names);

    /**
     * @brief Returns the set of aliases associated with this block.
     *
     * The returned set is a copy; modifications to it do not affect the
     * original BlockName instance.
     *
     * @return A copy of the internal alias set.
     */
    std::unordered_set<std::string> get_alias() const;

    /**
     * @brief Implicit conversion to std::string.
     *
     * If the block has exactly one alias, that alias is returned.
     * If multiple aliases are present, a warning is emitted via LOG_WARN
     * and one arbitrary alias (the first element of the internal set)
     * is returned, thus potentially discarding information.
     *
     * @return One alias representing this block (or an empty string if
     *         no aliases are stored).
     */
    operator std::string() const;

    /**
     * @brief Checks whether a given alias is registered for this block.
     *
     * @param alias Alias to look for.
     * @return true if @p alias is in the internal alias set, false otherwise.
     */
    bool hasAlias(const std::string& alias) const;

    /**
     * @brief Returns a string representation of the block name.
     *
     * Equivalent to the implicit conversion operator to std::string.
     * Mainly provided for explicit, non-implicit usage.
     *
     * @return One alias representing this block, or an empty string if none.
     */
    std::string to_string() const;

    /**
     * @brief Equality comparison between two BlockName objects.
     *
     * Two BlockName objects are considered equal if they share at least
     * one common alias. In particular, this means that:
     *   - equality is not simply based on the exact same alias set,
     *   - blocks with disjoint alias sets are considered different.
     *
     * @param other BlockName to compare with.
     * @return true if there exists at least one alias common to both blocks.
     */
    bool operator==(const BlockName& other) const;

    /**
     * @brief Inequality comparison between two BlockName objects.
     *
     * @param other BlockName to compare with.
     * @return true if the blocks do not share any common alias.
     */
    bool operator!=(const BlockName& other) const;

    /**
     * @brief Comparison between a BlockName and a std::string.
     *
     * Returns true if @p name is one of the aliases of this block.
     *
     * @param name Alias to check.
     * @return true if @p name is a known alias, false otherwise.
     */
    bool operator==(const std::string& name) const;

    /**
     * @brief Comparison between a BlockName and a C-style string.
     *
     * Equivalent to comparing with a std::string constructed from @p name.
     *
     * @param name Alias to check (C-string).
     * @return true if @p name is a known alias, false otherwise.
     */
    bool operator==(const char* name) const;

    /**
     * @brief Inequality comparison between a BlockName and a std::string.
     *
     * @param name Alias to check.
     * @return true if @p name is not a known alias of this block.
     */
    bool operator!=(const std::string& name) const;

    /**
     * @brief Adds a new alias to this block.
     *
     * Inserts @p alias into the internal alias set. If the alias already
     * exists, the set is unchanged.
     *
     * @param alias New alias to add.
     * @return Reference to the current BlockName (for chaining).
     */

    BlockName& addAlias(const std::string& alias);

    /**
     * @brief Concatenation operator between a prefix string and a BlockName.
     *
     * Creates a new BlockName whose aliases are obtained by prefixing
     * each alias of @p rhs with @p lhs.
     *
     * Example:
     * @code
     * BlockName b{"MASS","mass"};
     * BlockName c = "BLOCK_" + b; // aliases: {"BLOCK_MASS", "BLOCK_mass"}
     * @endcode
     *
     * @param lhs Prefix string.
     * @param rhs BlockName whose aliases will be prefixed.
     * @return A new BlockName with concatenated aliases.
     */
    friend BlockName operator+(const std::string& lhs, const BlockName& rhs);

    /**
     * @brief Stream output operator for BlockName.
     *
     * All aliases are written to the stream, separated by '/'.
     * The order is determined by the iteration order of the internal
     * unordered_set and is therefore not guaranteed to be stable.
     *
     * Example output: "SMINPUTS/sminputs".
     *
     * @param os Output stream.
     * @param block BlockName to print.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const BlockName& block);

    /**
     * @brief Converts all aliases to upper-case in-place.
     *
     * For each stored alias, all characters are transformed using std::toupper.
     * The resulting set of upper-case aliases replaces the original one.
     *
     * This is often used to enforce the conventional SLHA behavior that
     * block names are case-insensitive but typically stored in upper case.
     */
    void to_upper();

    /**
     * @brief Strict ordering operator for BlockName.
     *
     * The ordering is defined by:
     *   - first constructing two std::set<std::string> from the alias sets
     *     of @c *this and @p other, and
     *   - then performing the standard lexicographic comparison between
     *     those ordered sets.
     *
     * This makes BlockName usable as a key in ordered containers such as
     * std::map or std::set.
     *
     * @param other BlockName to compare with.
     * @return true if the sorted alias set of @c *this is lexicographically
     *         smaller than that of @p other.
     */
    bool operator<(const BlockName& other) const;
};

/**
 * @brief Symmetric comparison between std::string and BlockName.
 *
 * This overload allows code such as:
 * @code
 * if ("SMINPUTS" == block) { ... }
 * @endcode
 *
 * @param lhs Alias to compare.
 * @param rhs BlockName to compare against.
 * @return true if @p lhs is a known alias of @p rhs.
 */
inline bool operator==(const std::string& lhs, const BlockName& rhs) {
    return rhs == lhs;
}

/**
 * @brief Symmetric inequality between std::string and BlockName.
 *
 * @param lhs Alias to compare.
 * @param rhs BlockName to compare against.
 * @return true if @p lhs is not a known alias of @p rhs.
 */
inline bool operator!=(const std::string& lhs, const BlockName& rhs) {
    return !(rhs == lhs);
}

namespace std {
    /**
     * @brief Hash functor specialization for BlockName.
     *
     * This allows BlockName to be used as a key in std::unordered_map or
     * std::unordered_set. The hash is built from all aliases using a
     * Boost-style hash combine scheme in order to be sensitive to the
     * whole alias set (but not to its order).
     */
    template <>
    struct hash<BlockName> {
        /**
         * @brief Computes a hash value for a BlockName.
         *
         * Each alias contributes to the overall hash via a hash-combine-like
         * operation. Two BlockName objects that carry the same alias set
         * (ignoring ordering) will produce the same hash.
         *
         * @param p BlockName to hash.
         * @return Combined hash value.
         */
        std::size_t operator()(const BlockName& p) const noexcept {
            size_t h = 0;
            for (const auto& name : p.get_alias()) {
                h ^= std::hash<std::string>{}(name) + 0x9e3779b9 + (h << 6) + (h >> 2); // boost-style hash combine
            }
            return h;
        }
    };
}

#endif