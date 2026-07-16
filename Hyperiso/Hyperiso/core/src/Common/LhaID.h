#ifndef LHA_ID_H
#define LHA_ID_H

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
 * @file LhaID.h
 * @brief Lightweight representation of LHAPDF / SLHA-style identifiers.
 *
 * This header defines the LhaID struct, a small utility type used to encode
 * identifiers that may be composed of several integer components
 * (e.g. PDG-like IDs, block indices, etc.).
 *
 * The representation is essentially a vector of long integers, with helpers
 * to:
 *   - construct from integers, initializer lists, vectors or strings,
 *   - serialize to a compact string form (e.g. "1_2_3"),
 *   - compare and hash IDs for use in associative containers.
 */

/**
 * @struct LhaID
 * @brief Represents an identifier of a LHA element, possibly containing several sub-ids.
 *
 * An LhaID is internally stored as a sequence of integer parts. This allows
 * one to represent both simple one-dimensional identifiers (e.g. a PDG code)
 * and multi-dimensional indices (e.g. block/row/column) in a uniform way.
 *
 * The canonical string representation is obtained by joining the parts with
 * underscores, for example:
 *   - parts = { 511 }          → "511"
 *   - parts = { 5, 2 }         → "5_2"
 *   - parts = { 321, 1, 3 }    → "321_1_3"
 *
 * The type also supports:
 *   - comparison operators (lexicographic on parts),
 *   - a conversion to long for trivial (single-part) IDs, with a warning
 *     if information would be lost.
 */
struct LhaID {
    std::vector<long> parts;     /**< Collection of sub-ids. */

    /**
     * @brief Constructs an LhaID from a list of integral sub-ids.
     *
     * This variadic constructor accepts any number of arguments that are
     * convertible to long and stores them in order in the @ref parts vector.
     *
     * Examples:
     * @code
     * LhaID a(511);          // parts = {511}
     * LhaID b(5, 2);         // parts = {5, 2}
     * LhaID c(321, 1, 3);    // parts = {321, 1, 3}
     * @endcode
     *
     * @tparam Args Parameter pack of types convertible to long.
     * @param sub_ids Sub-ids of the element.
     */
    template<typename... Args>
    requires (std::convertible_to<Args, long> && ...)
    LhaID(Args... sub_ids) : parts({static_cast<long>(sub_ids)...}) {}

    /**
     * @brief Constructs an LhaID from a string of underscore-separated integers.
     *
     * The string is split on the '_' character and each token is converted
     * to a long using std::stol. An empty string corresponds to an empty
     * parts vector.
     *
     * Examples:
     *  - "511"      → parts = {511}
     *  - "5_2"      → parts = {5, 2}
     *  - ""         → parts = {}
     *
     * @param parts String encoding the sub-ids of the element, separated by '_'.
     *
     * @throws std::invalid_argument or std::out_of_range if a token cannot
     *         be converted by std::stol.
     */
    LhaID(const std::string& parts);

    /**
     * @brief Constructs an LhaID from an explicit vector of sub-ids.
     *
     * The content of @p sub_ids is moved into the internal @ref parts vector.
     *
     * @param sub_ids Sub-identifiers of the element.
     */
    LhaID(const std::vector<long>& sub_ids) : parts(sub_ids) {}
    LhaID(std::vector<long>&& sub_ids) : parts(std::move(sub_ids)) {}

    // /**
    //  * @brief Constructs an LhaID from an initializer list of sub-ids.
    //  *
    //  * This is especially convenient for brace-initialization:
    //  * @code
    //  * LhaID id{5, 2};
    //  * @endcode
    //  *
    //  * @param sub_ids Sub-identifiers of the element.
    //  */
    // LhaID(std::initializer_list<long> sub_ids) : parts(sub_ids) {}

    /**
     * @brief Constructs an LhaID with a single identifier.
     *
     * Equivalent to calling the variadic constructor with a single argument
     * or using an initializer list with one element.
     *
     * @param id Identifier of the element.
     */
    LhaID(long id) : parts({id}) {}

    /**
     * @brief Returns the canonical string representation of this LhaID.
     *
     * The parts are joined with underscores in order. If @ref parts is empty,
     * an empty string is returned.
     *
     * @return String representation, e.g. "5_2".
     */
    std::string to_string() const;

    /**
     * @brief Returns the underlying vector of sub-ids.
     *
     * This is a value-returning accessor; modifications of the returned
     * vector do not affect the original LhaID.
     *
     * @return A copy of the internal parts vector.
     */
    std::vector<long> get_parts() const { return parts; }
    
    /**
     * @brief Implicit conversion of a trivial LhaID to a long integer.
     *
     * This conversion is intended for IDs that contain a single part
     * (e.g. simple PDG codes). If @ref parts contains more than one element,
     * a warning is emitted via LOG_WARN and only the first element is
     * returned, effectively discarding information.
     *
     * @return The first element of @ref parts.
     *
     * @warning For multi-part IDs, only the first component is kept.
     */
    operator long() const;

    /**
     * @brief Equality comparison operator.
     *
     * Two LhaID objects are equal if and only if their @ref parts vectors
     * are equal.
     *
     * @param lhs Left-hand side.
     * @param rhs Right-hand side.
     * @return true if both identifiers have identical parts, false otherwise.
     */
    inline friend bool operator==(const LhaID& lhs, const LhaID& rhs) { return lhs.parts == rhs.parts; };

    /**
     * @brief Inequality comparison operator.
     *
     * @param lhs Left-hand side.
     * @param rhs Right-hand side.
     * @return true if @p lhs and @p rhs differ, false otherwise.
     */
    inline friend bool operator!=(const LhaID& lhs, const LhaID& rhs) { return !(lhs == rhs); };

    /**
     * @brief Strict weak ordering based on lexicographic comparison of parts.
     *
     * This operator enables LhaID to be used as a key in ordered containers
     * such as std::map or std::set.
     *
     * @param lhs Left-hand side.
     * @param rhs Right-hand side.
     * @return true if @p lhs is lexicographically smaller than @p rhs.
     */
    inline friend bool operator<(const LhaID& lhs, const LhaID& rhs) { return (lhs.parts <=> rhs.parts) == std::weak_ordering::less; };

    /**
     * @brief Stream output operator for LhaID.
     *
     * Writes the canonical string representation (as returned by
     * to_string()) to the output stream.
     *
     * @param os Output stream.
     * @param id Identifier to print.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream&, const LhaID&);
};

namespace std {
    /**
     * @brief Hash functor specialization for LhaID.
     *
     * This specialization allows LhaID to be used as a key in std::unordered_map
     * or std::unordered_set. The hash is computed from the canonical string
     * representation of the identifier.
     */
    template <>
    struct hash<LhaID> {
        /**
         * @brief Computes the hash of a given LhaID.
         *
         * @param p LhaID to hash.
         * @return Hash value derived from p.to_string().
         */
        std::size_t operator()(const LhaID& p) const noexcept {
            return std::hash<std::string>{}(p.to_string());
        }
    };
}

#endif