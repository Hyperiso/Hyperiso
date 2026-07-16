#ifndef ENUM_MAPPER_H
#define ENUM_MAPPER_H

#include <map>
#include <string>
#include <vector>

/**
 * @file EnumMapper.h
 * @brief Generic utilities for bidirectional enum–string mappings.
 *
 * This header provides:
 *   - a CRTP base class EnumMapperBase for defining enum-to-string and
 *     string-to-enum mappers,
 *   - a small helper function invert_map() to build inverse maps.
 *
 * A typical usage pattern is:
 *   - define an enum type,
 *   - define a struct Derived : EnumMapperBase<Enum, Derived> that exposes:
 *       * static const std::map<Enum, std::string>& mapping();
 *       * static const std::map<std::string, Enum>& inverse_mapping();
 *   - use Derived::str(e) to get the name of an enum value,
 *     and Derived::enum_elt(name) to recover the enum from its string.
 */

/**
 * @class EnumMapperBase
 * @brief CRTP base class for enum–string mappers.
 *
 * This template implements a small set of convenience functions relying on
 * a Derived class that must provide:
 *   - static const std::map<Enum, std::string>& mapping();
 *   - static const std::map<std::string, Enum>& inverse_mapping();
 *
 * The design is based on CRTP (Curiously Recurring Template Pattern):
 *   struct MyEnumMap : EnumMapperBase<MyEnum, MyEnumMap> { ... };
 *
 * The implementation assumes that:
 *   - mapping() returns an ordered map that is stable across calls,
 *   - inverse_mapping() is consistent with mapping(), typically built via
 *     invert_map(mapping()).
 *
 * @tparam Enum    Enum type being mapped.
 * @tparam Derived CRTP derived type implementing mapping() / inverse_mapping().
 */
template <typename Enum, typename Derived>
class EnumMapperBase {
public:
    /**
     * @brief Returns the string representation of an enum value.
     *
     * This function simply forwards to Derived::mapping().at(e). If @p e
     * is not present in the mapping, std::out_of_range is thrown.
     *
     * @param e Enum value to convert.
     * @return Corresponding string name.
     *
     * @throws std::out_of_range if the enum value is not in the mapping.
     */
    static std::string str(Enum e) {
        return Derived::mapping().at(e);
    }

    /**
     * @brief Returns the enum value associated to a given string.
     *
     * This function forwards to Derived::inverse_mapping().at(name).
     * If @p name is not present in the inverse mapping, std::out_of_range
     * is thrown.
     *
     * @param name String representation of the enum value.
     * @return Corresponding enum value.
     *
     * @throws std::out_of_range if the string is not in the inverse mapping.
     */
    static Enum enum_elt(const std::string& name) {
        return Derived::inverse_mapping().at(name);
    }

    /**
     * @brief Returns all string representations defined in the mapping.
     *
     * The strings are returned in the iteration order of Derived::mapping(),
     * which is the key order of the underlying std::map<Enum, std::string>.
     *
     * @return Vector containing all mapped strings.
     */
    static std::vector<std::string> get_str() {
        std::vector<std::string> out;
        for (const auto& [key, val] : Derived::mapping()) {
            out.push_back(val);
        }
        return out;
    }

    /**
     * @brief Returns all enum values defined in the mapping.
     *
     * The enum values are returned in the iteration order of Derived::mapping(),
     * i.e. the key order of the underlying std::map<Enum, std::string>.
     *
     * @return Vector containing all mapped enum values.
     */
    static std::vector<Enum> get_enum() {
        std::vector<Enum> out;
        for (const auto& [key, val] : Derived::mapping()) {
            out.push_back(key);
        }
        return out;
    }
};

/**
 * @brief Builds the inverse of a std::map by swapping keys and values.
 *
 * Given an ordered map @p original of type std::map<K, V>, this helper
 * constructs a new std::map<V, K> where each value becomes a key and each
 * key becomes a value.
 *
 * If multiple keys in @p original share the same value, the last one
 * encountered in iteration order will overwrite the previous entry in
 * the inverse map. This is intentional and tested behaviour.
 *
 * @tparam K Key type of the original map.
 * @tparam V Value type of the original map.
 *
 * @param original Map to invert.
 * @return Inverted map where values become keys and keys become values.
 */
template <typename K, typename V>
std::map<V, K> invert_map(const std::map<K, V>& original) {
    std::map<V, K> inv;
    for (const auto& [k, v] : original) {
        inv[v] = k;
    }
    return inv;
}

#endif