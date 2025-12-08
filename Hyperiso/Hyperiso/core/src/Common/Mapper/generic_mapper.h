#ifndef GENERIC_MAPPER_H
#define GENERIC_MAPPER_H

#include <map>
#include <optional>
#include <string_view>
#include <vector>
#include "dynamic_registry.h"
#include "Map.h"

/**
 * @file generic_mapper.h
 * @brief Generic enum↔string (and optional external key) mappers built on DynamicRegistry.
 *
 * This header provides two templated mapper classes:
 *   - GenericMapperNoExt:  bridges an enum type and its string names, with
 *                          support for dynamic aliases and user-defined symbols.
 *   - GenericMapperWithExt: extends the above with an additional "external key"
 *                           (e.g. FLHA indices, integer codes, etc).
 *
 * Both mappers are parameterized by:
 *   - Tag    : a tag type used to instantiate SymbolId<Tag> (strong-typed string id),
 *   - EnumT  : the enum type being mapped,
 *   - MapFn  : a function pointer returning a const std::map<EnumT,std::string>&
 *              with the builtin mapping (legacy/static map),
 *   - ExtFn  : (for GenericMapperWithExt) a function pointer returning
 *              const std::map<EnumT,ExternalKey>& for the external mapping.
 *
 * The internal DynamicRegistry uses case-insensitive normalization via
 * normalize_key(), so lookups are case-insensitive and aliases can be
 * added at runtime without interfering with builtin names.
 */

// Utility
/**
 * @brief Convenience alias for the identifier type associated with a given tag.
 *
 * @tparam Tag Tag used for SymbolId instantiation.
 */
template<class Tag>
using IdOf = SymbolId<Tag>;

/**
 * @brief Shortcut helper to normalize a string_view into a canonical key.
 *
 * This simply forwards to normalize_key().
 *
 * @param s Input string_view.
 * @return Normalized string.
 */
inline std::string nk(std::string_view s){ return normalize_key(s); }

// ------------------------------ Without External Keys ------------------------------

/**
 * @class GenericMapperNoExt
 * @brief Generic enum↔string mapper without external keys.
 *
 * This template provides a high-level interface to:
 *   - expose builtin mappings between EnumT and std::string (via MapFn),
 *   - support case-insensitive lookups,
 *   - register additional user-defined symbols and aliases at runtime,
 *   - map back and forth between EnumT and SymbolId<Tag>.
 *
 * Initialization:
 *   - Builtins are registered lazily on first use via ensure_init(), which
 *     reads the map returned by MapFn() and inserts entries into an internal
 *     DynamicRegistry<Tag, void, void>.
 *   - init_builtins() can be called explicitly to force initialization.
 *
 * Typical usage pattern:
 * @code
 * // Suppose observable_mapping() returns const std::map<Observables,std::string>&
 * struct ObsTag {};
 * using ObsMapper = GenericMapperNoExt<ObsTag, Observables, observable_mapping>;
 *
 * auto names = ObsMapper::get_str();            // builtin names
 * auto all   = ObsMapper::get_str_all();        // builtin + custom names
 * auto id    = ObsMapper::id_of("BR(B->Kll)");  // lookup by name (case-insensitive)
 * auto s     = ObsMapper::str(Observables::SomeObs);
 * @endcode
 *
 * @tparam Tag    Tag used for the associated SymbolId<Tag>.
 * @tparam EnumT  Enum type being mapped.
 * @tparam MapFn  Function returning the builtin map EnumT→std::string.
 */
template<class Tag, class EnumT,
         const std::map<EnumT, std::string>& (*MapFn)()>
class GenericMapperNoExt {
    using Reg = DynamicRegistry<Tag, void, void>;

    /// Returns the singleton registry instance.
    static Reg& reg();

    /// Returns a reference to the static "initialized" flag.
    static bool& inited();

    /// Ensures builtin entries are registered in the registry (lazy init).
    static void ensure_init();
public:
    /**
     * @brief Explicitly initializes builtin entries.
     *
     * This simply forwards to ensure_init(). It can be used at startup
     * to populate the registry eagerly instead of lazily.
     */
    static void init_builtins();

    /**
     * @brief Legacy enum lookup using only the static MapFn() table.
     *
     * This function emulates the pre-registry behavior by scanning the
     * map returned by MapFn() and matching against the normalized key.
     * It does not consult the dynamic registry and therefore ignores
     * custom registrations.
     *
     * @param s Name to look up (case-insensitive).
     * @return Matching EnumT value.
     *
     * @throws std::out_of_range if no matching enum value is found.
     */
    static EnumT enum_elt_legacy(std::string_view s);
    
    /**
     * @brief Returns the list of builtin string names from MapFn().
     *
     * Only entries in the static map are returned; custom/user-defined
     * aliases are not included.
     *
     * @return Vector of builtin names.
     */
    static std::vector<std::string> get_str();

    /**
     * @brief Returns the list of builtin enum values from MapFn().
     *
     * @return Vector of EnumT values as defined in the static map.
     */
    static std::vector<EnumT> get_enum();

    /**
     * @brief Returns all known names: builtin and custom.
     *
     * This triggers initialization and then queries the underlying
     * DynamicRegistry for all registered identifiers.
     *
     * @return Vector of all distinct identifier strings.
     */
    static std::vector<std::string> get_str_all();

    /**
     * @brief Resolves a string into an IdOf<Tag> via the registry.
     *
     * The lookup is case-insensitive and includes both builtin and
     * custom-registered entries.
     *
     * @param s Name or alias to look up.
     * @return Corresponding IdOf<Tag>.
     *
     * @throws std::runtime_error if the name is unknown to the registry.
     */
    static IdOf<Tag> id_of(std::string_view s);

    /**
     * @brief Alias for id_of(), kept for compatibility.
     *
     * @param s Name or alias to look up.
     * @return Corresponding IdOf<Tag>.
     */
    static IdOf<Tag> enum_elt(std::string_view s);

    /**
     * @brief Registers a custom (non-builtin) entry in the registry.
     *
     * The canonical name is stored as the primary identifier, and the
     * supplied aliases are also registered as case-insensitive keys
     * mapping to the same IdOf<Tag>.
     *
     * @param canonical Canonical name for the new identifier.
     * @param aliases   Optional list of aliases for lookup.
     * @return true if registration succeeded, false if blocked by a
     *         builtin entry with the same canonical name.
     */
    static bool register_custom(std::string canonical, std::vector<std::string> aliases={});

    /**
     * @brief Lists all distinct IdOf<Tag> entries in the registry.
     *
     * This includes both builtin and custom entries, each appearing
     * only once even if it has multiple aliases.
     *
     * @return Vector of IdOf<Tag>.
     */
    static std::vector<IdOf<Tag>> list_all();

    /**
     * @brief Returns the string representation associated with an identifier.
     *
     * This simply returns id.str(), i.e. the underlying canonical name
     * stored in the SymbolId<Tag>.
     *
     * @param id Identifier.
     * @return String representation of @p id.
     */
    static std::string str(const IdOf<Tag>& id);

    /**
     * @brief Returns the builtin string name corresponding to an enum value.
     *
     * @param e Enum value.
     * @return Corresponding string name from MapFn().
     *
     * @throws std::out_of_range if @p e is not present in the static map.
     */
    static std::string str(EnumT e);

    /**
     * @brief Converts an enum value to an IdOf<Tag>.
     *
     * The resulting identifier wraps the builtin string name obtained
     * from MapFn().at(e).
     *
     * @param e Enum value.
     * @return Corresponding IdOf<Tag>.
     */
    static IdOf<Tag> to_id(EnumT e);

    /**
     * @brief Attempts to recover the enum value associated with an identifier.
     *
     * The identifier's string is normalized and compared against the
     * builtin names from MapFn(). This does not consult the dynamic
     * registry, so only builtin names are considered.
     *
     * @param id Identifier to decode.
     * @return Matching EnumT if found, std::nullopt otherwise.
     */
    static std::optional<EnumT> enum_of(const IdOf<Tag>& id);
};

// ------------------------------ With External Keys ------------------------------

/**
 * @class GenericMapperWithExt
 * @brief Generic enum↔string mapper with external keys and dynamic aliases.
 *
 * This template extends GenericMapperNoExt by adding support for
 * an "external key" (ExternalKey), which can be used to represent
 * e.g. FLHA indices, integer codes, or other external identifiers.
 *
 * Initialization:
 *   - Builtin entries are created from the map returned by MapFn() and
 *     the external-key map returned by ExtFn().
 *   - Each EnumT value is mapped to a canonical name and optionally
 *     to an external key, stored in a DynamicRegistry<Tag,ExternalKey,Hash>.
 *
 * Features:
 *   - Bidirectional mapping between EnumT and IdOf<Tag>.
 *   - Bidirectional mapping between external keys and IdOf<Tag>.
 *   - Case-insensitive lookup and dynamic alias registration.
 *
 * @tparam Tag         Tag used for SymbolId<Tag>.
 * @tparam EnumT       Enum type to map.
 * @tparam ExternalKey Type for the external key.
 * @tparam Hash        Hash functor for ExternalKey in unordered_map.
 * @tparam MapFn       Function returning builtin EnumT→std::string map.
 * @tparam ExtFn       Function returning EnumT→ExternalKey map.
 */
template<class Tag, class EnumT, class ExternalKey, class Hash,
         const std::map<EnumT, std::string>& (*MapFn)(),
         const std::map<EnumT, ExternalKey>& (*ExtFn)()>
class GenericMapperWithExt {
    using Reg = DynamicRegistry<Tag, ExternalKey, Hash>;

    /// Returns the singleton registry instance.
    static Reg& reg();

    /// Returns a reference to the static "initialized" flag.
    static bool& inited();

    /// Ensures builtin entries are registered in the registry (lazy init).
    static void ensure_init();
public:
    /**
     * @brief Explicitly initializes builtin entries.
     *
     * Populates the registry from the static maps provided by MapFn()
     * and ExtFn() if not already initialized.
     */
    static void init_builtins();

    /**
     * @brief Returns the list of builtin string names from MapFn().
     *
     * Only the static mapping is consulted; dynamic/custom entries
     * are not included.
     *
     * @return Vector of builtin names.
     */
    static std::vector<std::string> get_str();

    /**
     * @brief Returns the list of builtin enum values.
     *
     * @return Vector of EnumT values from MapFn().
     */
    static std::vector<EnumT> get_enum();

    /**
     * @brief Returns all known names (builtin + custom).
     *
     * This triggers initialization and then enumerates the registry
     * to collect all identifier strings.
     *
     * @return Vector of all identifier names.
     */
    static std::vector<std::string> get_str_all();

    /**
     * @brief Legacy enum lookup using only the static MapFn() table.
     *
     * As in GenericMapperNoExt::enum_elt_legacy, this scans the static
     * map and ignores custom registrations and external keys.
     *
     * @param s Name to look up (case-insensitive).
     * @return Matching EnumT value.
     *
     * @throws std::out_of_range if no match is found.
     */
    static EnumT enum_elt_legacy(std::string_view s);


    /**
     * @brief Resolves a string into an IdOf<Tag> via the registry.
     *
     * The lookup is case-insensitive and takes into account builtin
     * and custom aliases.
     *
     * @param s Name or alias.
     * @return Corresponding IdOf<Tag>.
     *
     * @throws std::runtime_error if the name is unknown.
     */
    static IdOf<Tag> id_of(std::string_view s);

    /**
     * @brief Alias for id_of(), kept for compatibility.
     *
     * @param s Name or alias.
     * @return Corresponding IdOf<Tag>.
     */
    static IdOf<Tag> enum_elt(std::string_view s);

    /**
     * @brief Resolves an external key into an identifier.
     *
     * Uses the registry's external-key index to find the canonical
     * identifier associated with @p k.
     *
     * @param k External key to look up.
     * @return Matching IdOf<Tag>, or std::nullopt if none.
     */
    static std::optional<IdOf<Tag>> from_external(const ExternalKey& k);

    /**
     * @brief Returns the external key associated with a given identifier.
     *
     * @param id Identifier to query.
     * @return Matching ExternalKey if present, std::nullopt otherwise.
     */
    static std::optional<ExternalKey> external_of(const IdOf<Tag>& id);

    /**
     * @brief Registers a custom entry with optional external key.
     *
     * Similar to register_custom() in GenericMapperNoExt, but also
     * allows attaching an external key that participates in
     * from_external()/external_of() mappings.
     *
     * @param canonical Canonical name for the identifier.
     * @param aliases   Optional aliases for lookup.
     * @param ext       Optional external key to associate.
     * @return true if registration succeeded, false if blocked by
     *         an existing builtin entry.
     */
    static bool register_custom(std::string canonical, std::vector<std::string> aliases={}, std::optional<ExternalKey> ext=std::nullopt);

    /**
     * @brief Sets or updates the external key of an existing identifier.
     *
     * If the identifier exists in the registry, its external key is set
     * to @p k and the external-key index updated accordingly.
     *
     * @param id Identifier whose external key should be set.
     * @param k  New external key value.
     * @return true if the identifier was found and updated, false otherwise.
     */
    static bool set_external(const IdOf<Tag>& id, const ExternalKey& k);

    /**
     * @brief Lists all distinct identifiers currently registered.
     *
     * Includes both builtin and custom entries; each canonical id is
     * returned once, regardless of the number of aliases.
     *
     * @return Vector of IdOf<Tag>.
     */
    static std::vector<IdOf<Tag>> list_all();

    /**
     * @brief Returns the string representation of an identifier.
     *
     * This simply returns id.str().
     *
     * @param id Identifier.
     * @return String representation of @p id.
     */
    static std::string str(const IdOf<Tag>& id);

    /**
     * @brief Returns the builtin string name of an enum value.
     *
     * @param e Enum value.
     * @return Corresponding string from MapFn().
     *
     * @throws std::out_of_range if @p e is not present in the static map.
     */
    static std::string str(EnumT e);

    /**
     * @brief Converts an enum value to an IdOf<Tag>.
     *
     * The identifier wraps the builtin string returned by MapFn().at(e).
     *
     * @param e Enum value.
     * @return Corresponding IdOf<Tag>.
     */
    static IdOf<Tag> to_id(EnumT e);

    /**
     * @brief Attempts to recover the enum value associated with an identifier.
     *
     * The identifier’s string is normalized and compared against the
     * builtin names in MapFn(). This does not consult the dynamic
     * registry, so only builtin names are considered.
     *
     * @param id Identifier to decode.
     * @return Matching EnumT if found, std::nullopt otherwise.
     */
    static std::optional<EnumT> enum_of(const IdOf<Tag>& id);
};

#include "generic_mapper.tpp"

#endif