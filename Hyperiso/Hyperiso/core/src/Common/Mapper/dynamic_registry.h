#ifndef DYNAMIC_REGISTRY_H
#define DYNAMIC_REGISTRY_H

#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <optional>
#include <algorithm>
#include <cctype>

/**
 * @file SymbolRegistry.h
 * @brief Strong-typed string identifiers and case-insensitive dynamic registries.
 *
 * This header defines:
 *   - SymbolId<Tag>: a small wrapper around std::string providing strong typing
 *     via a tag type, with case-insensitive comparison and hashing.
 *   - DynamicRegistry<Tag, ExternalKey, Hash>: a generic registry used to
 *     associate SymbolId values with optional external keys and aliases.
 *
 * The registry is case-insensitive: all lookup keys and identifiers are
 * normalized via normalize_key(), which lowercases characters (and can
 * optionally ignore punctuation if desired).
 */

// ------------------------------------------------------------------
// SymbolId<Tag> : wrapper string-strong-typed
// ------------------------------------------------------------------

/**
 * @class SymbolId
 * @brief Strongly-typed string identifier.
 *
 * SymbolId<Tag> wraps a std::string in a lightweight type that is
 * distinguished by its tag parameter. This avoids mixing identifiers
 * that share the same underlying representation but belong to different
 * domains (e.g. observable IDs vs. parameter IDs).
 *
 * The class itself only stores the original string and exposes it via
 * str(). Equality, ordering and hashing are defined externally and
 * are case-insensitive via normalize_key().
 *
 * @tparam Tag Empty tag type used to distinguish different identifier kinds.
 */
template<class Tag>
class SymbolId {
    std::string name_;
public:
    /**
     * @brief Default constructor, creates an empty identifier.
     */
    SymbolId() = default;

    /**
     * @brief Constructs a SymbolId from a string.
     *
     * The string is stored as-is; normalization (lowercasing, etc.) is
     * performed only in comparison and hashing functions.
     *
     * @param n Identifier name.
     */
    explicit SymbolId(std::string n) : name_(std::move(n)) {}

    /**
     * @brief Returns the underlying string.
     *
     * This is the original, non-normalized name.
     *
     * @return Const reference to the stored string.
     */
    const std::string& str() const { return name_; }
};





// ------------------------------------------------------------------
// case-insensitive
// ------------------------------------------------------------------

/**
 * @brief Normalizes a string key in a case-insensitive manner.
 *
 * The current implementation:
 *   - lowercases all characters using std::tolower,
 *   - preserves all characters (alphanumeric and others).
 *
 * Optionally, the commented-out line can be enabled to ignore characters
 * that are not alphanumeric (such as spaces or punctuation).
 *
 * @param s Input string view to normalize.
 * @return Lowercased std::string, suitable for use as a canonical key.
 */
inline std::string normalize_key(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s) {
        // if (!std::isalnum(c)) continue; // <-- option : ignore space and punctuation
        out.push_back(static_cast<char>(std::tolower(c)));
    }
    return out;
}

/**
 * @brief Case-insensitive equality comparison for SymbolId<Tag>.
 *
 * Two SymbolId<Tag> values are considered equal if their normalized
 * names (as returned by normalize_key()) are equal, independently of
 * original case.
 *
 * @tparam Tag Tag type of the SymbolId.
 * @param a First identifier.
 * @param b Second identifier.
 * @return true if both identifiers are equal under normalization.
 */
template<class Tag>
inline bool operator==(const SymbolId<Tag>& a, const SymbolId<Tag>& b) noexcept {
    return normalize_key(a.str()) == normalize_key(b.str());
}

/**
 * @brief Case-insensitive strict ordering for SymbolId<Tag>.
 *
 * The ordering is defined by lexicographic comparison of the normalized
 * names (see normalize_key()). This allows SymbolId<Tag> to be used as
 * keys in ordered containers if needed.
 *
 * @tparam Tag Tag type of the SymbolId.
 * @param a First identifier.
 * @param b Second identifier.
 * @return true if the normalized name of @p a is lexicographically
 *         smaller than that of @p b.
 */
template<class Tag>
inline bool operator<(const SymbolId<Tag>& a, const SymbolId<Tag>& b) noexcept {
    const auto an = normalize_key(a.str());
    const auto bn = normalize_key(b.str());
    return an < bn;
}

 /**
 * @brief Hash specialization for SymbolId<Tag>.
 *
 * SymbolId<Tag> is hashed by first normalizing its string via
 * normalize_key() and then using std::hash<std::string>. This makes
 * the hash consistent with the case-insensitive equality operator.
 *
 * @tparam Tag Tag type of the SymbolId.
 */
namespace std {
    template<class Tag>
    struct hash<SymbolId<Tag>> {
        /**
         * @brief Computes the hash value of a SymbolId<Tag>.
         *
         * @param id Identifier to hash.
         * @return Hash of the normalized string.
         */
        size_t operator()(const SymbolId<Tag>& id) const noexcept {
            return std::hash<std::string>{}(normalize_key(id.str()));
        }
    };
}


// ------------------------------------------------------------------
// DynamicRegistry: registry of SymbolId with optional external key
// ------------------------------------------------------------------

/**
 * @class DynamicRegistry
 * @brief Case-insensitive registry of SymbolId entries with optional external keys.
 *
 * The primary template manages:
 *   - a mapping from normalized string names (and aliases) to Entry objects,
 *   - an optional mapping from external keys (e.g. integer codes) to names.
 *
 * Typical usage:
 * @code
 * struct MyTag {};
 * using MyId = SymbolId<MyTag>;
 * using MyRegistry = DynamicRegistry<MyTag, int>;
 *
 * MyRegistry reg;
 * reg.register_id(MyId{"Foo"}, {"foo", "FOO"}, 42, true);
 *
 * auto id = reg.find("fOo"); // returns MyId{"Foo"}
 * auto id2 = reg.find_by_external(42);
 * @endcode
 *
 * Builtin entries can be protected against being overridden by later
 * registrations.
 *
 * @tparam Tag         Tag used for the SymbolId type.
 * @tparam ExternalKey Type of the external key (e.g. int, std::string).
 * @tparam Hash        Hash functor type used for ExternalKey.
 */
template<class Tag, class ExternalKey, class Hash = std::hash<ExternalKey>>
class DynamicRegistry {
public:
    /**
     * @brief Entry stored in the registry.
     *
     * Each Entry contains:
     *   - an identifier (id),
     *   - a list of aliases (all case-insensitive),
     *   - an optional external key,
     *   - a flag indicating whether the entry is builtin.
     */
    struct Entry {
        SymbolId<Tag> id;                   ///< Canonical identifier.
        std::vector<std::string> aliases;   ///< Additional names mapping to this entry.
        std::optional<ExternalKey> ext;     ///< Optional external key.
        bool builtin = false;               ///< If true, entry cannot be overridden
    };

private:
    /// Map from normalized name to Entry.
    std::unordered_map<std::string, Entry> by_name_;
    /// Map from external key to normalized canonical name.
    std::unordered_map<ExternalKey, std::string, Hash> by_ext_;

public:
    /**
     * @brief Registers or updates an identifier in the registry.
     *
     * The given @p id is inserted with the provided aliases and optional
     * external key. The canonical map key is normalize_key(id.str()).
     *
     * Behavior:
     *   - If the canonical name does not exist yet, a new entry is created.
     *   - If it exists and is marked as builtin, the registration fails
     *     and the function returns false.
     *   - Otherwise, the existing entry is replaced.
     *   - All aliases are also indexed (case-insensitive), pointing to
     *     the same Entry, unless they belong to a builtin entry that
     *     should not be overwritten.
     *   - If an external key is provided, it is registered in by_ext_.
     *
     * @param id      Identifier to register.
     * @param aliases List of aliases for lookup.
     * @param ext     Optional external key to associate.
     * @param builtin Whether this is a builtin entry (protected).
     * @return true if the registration succeeded, false if blocked by
     *         an existing builtin entry.
     */
    bool register_id(const SymbolId<Tag>& id,
                     std::vector<std::string> aliases,
                     std::optional<ExternalKey> ext,
                     bool builtin);

    /**
     * @brief Finds an identifier by name (or alias), case-insensitively.
     *
     * The input @p name is normalized via normalize_key() and looked up
     * in the internal map. If found, the corresponding SymbolId<Tag>
     * is returned.
     *
     * @param name Name or alias to search for.
     * @return The corresponding SymbolId<Tag>, or std::nullopt if not found.
     */
    std::optional<SymbolId<Tag>> find(std::string_view name) const;

    /**
     * @brief Finds an identifier by external key.
     *
     * Looks up @p k in the external-key map and then resolves the associated
     * canonical name to retrieve the SymbolId<Tag>.
     *
     * @param k External key value.
     * @return The corresponding SymbolId<Tag>, or std::nullopt if not found.
     */
    std::optional<SymbolId<Tag>> find_by_external(const ExternalKey& k) const;

    /**
     * @brief Returns the external key associated with a given identifier.
     *
     * @param id SymbolId to query.
     * @return The associated external key, or std::nullopt if none.
     */
    std::optional<ExternalKey> external_of(const SymbolId<Tag>& id) const;

    /**
     * @brief Updates or sets the external key of a given identifier.
     *
     * If the identifier is found in the registry, its external key is
     * set to @p k and the external-key index is updated accordingly.
     *
     * @param id Identifier whose external key is to be set.
     * @param k  New external key value.
     * @return true if the identifier exists and was updated, false otherwise.
     */
    bool update_external(const SymbolId<Tag>& id, const ExternalKey& k);

    /**
     * @brief Lists all distinct identifiers stored in the registry.
     *
     * The registry may contain multiple entries for different aliases
     * pointing to the same canonical id. This function returns each
     * distinct canonical id exactly once.
     *
     * @return Vector of all unique SymbolId<Tag> values.
     */
    std::vector<SymbolId<Tag>> list_all() const;
};

// ------------------------------------------------------------------
// DynamicRegistry specialization without external key
// ------------------------------------------------------------------

/**
 * @class DynamicRegistry<Tag, void, void>
 * @brief Case-insensitive registry of SymbolId entries without external keys.
 *
 * This partial specialization handles the case where no external key
 * is needed. It maintains only a name/alias-based registry.
 *
 * Usage:
 * @code
 * struct ObsTag{};
 * using ObsId = SymbolId<ObsTag>;
 * using ObsRegistry = DynamicRegistry<ObsTag, void, void>;
 *
 * ObsRegistry reg;
 * reg.register_id(ObsId{"BR_B2K"}, {"b->k", "B2K"}, true);
 * auto id = reg.find("b2k");
 * @endcode
 *
 * @tparam Tag Tag used for the SymbolId type.
 */
template<class Tag>
class DynamicRegistry<Tag, void, void> {
public:
    /**
     * @brief Entry stored in the registry (no external key).
     */
    struct Entry {
        SymbolId<Tag> id;                   ///< Canonical identifier.
        std::vector<std::string> aliases;   ///< Alias names.
        bool builtin = false;               ///< If true, protected from replacement.
    };
private:
    /// Map from normalized name/alias to Entry.
    std::unordered_map<std::string, Entry> by_name_;

public:
    /**
     * @brief Registers or updates an identifier (no external key).
     *
     * Semantics are similar to the primary template but without any external
     * key handling:
     *   - a new Entry is inserted or replaces a non-builtin one,
     *   - builtin entries cannot be replaced,
     *   - all aliases are indexed as case-insensitive keys.
     *
     * @param id      Identifier to register.
     * @param aliases List of aliases for lookup.
     * @param builtin Whether this is a builtin entry (protected).
     * @return true if registration succeeded, false if blocked by a builtin.
     */
    bool register_id(const SymbolId<Tag>& id,
                     std::vector<std::string> aliases,
                     bool builtin);

    /**
     * @brief Finds an identifier by name (or alias), case-insensitively.
     *
     * @param name Name or alias to search for.
     * @return The corresponding SymbolId<Tag>, or std::nullopt if not found.
     */
    std::optional<SymbolId<Tag>> find(std::string_view name) const;

    /**
     * @brief Lists all distinct identifiers in the registry.
     *
     * Duplicate entries arising from aliases are filtered out so that each
     * canonical id appears exactly once in the returned vector.
     *
     * @return Vector of unique SymbolId<Tag> values.
     */
    std::vector<SymbolId<Tag>> list_all() const;
};

#include "dynamic_registry.tpp"

#endif