#ifndef HYPERISO_UTILS_H
#define HYPERISO_UTILS_H

#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <memory>
#include <limits>

/**
 * @file HyperIsoUtils.h
 * @brief Miscellaneous utility types and functions for numerical and container handling.
 *
 * This header gathers a set of small, reusable utilities that are used
 * throughout the codebase:
 *   - aliases and constants for numerical work,
 *   - small string helpers (suffix check, splitting, lowercase conversion),
 *   - helpers for extracting map keys and comparing key sets,
 *   - generic printing utilities for values and maps,
 *   - simple caching and linear interpolation helpers for functions
 *     on a finite interval.
 */

 /**
 * @typedef complex_t
 * @brief Convenience alias for std::complex<double>.
 */
typedef std::complex<double> complex_t;

/**
 * @brief Numeric limits for double-precision floating point.
 *
 * This object is used to access properties such as machine epsilon
 * via nld.epsilon() in numerical routines (e.g. cache filling).
 */
const std::numeric_limits<double> nld = *new std::numeric_limits<double>;

/**
 * @brief Tests whether a string ends with a given suffix.
 *
 * The comparison is case-sensitive and operates on the raw bytes of the
 * strings.
 *
 * @param str    String to test.
 * @param suffix Suffix to look for.
 * @return true if @p str ends with @p suffix, false otherwise.
 */
bool ends_with(const std::string& str, const std::string& suffix);

/**
 * @brief Splits a string into parts using a single-character delimiter.
 *
 * Consecutive delimiters yield empty tokens. The delimiter itself is
 * not included in the resulting strings.
 *
 * Example:
 *   - split("a,b,c", ',')   → {"a","b","c"}
 *   - split("a,,c", ',')    → {"a","","c"}
 *
 * @param s         Input string to split.
 * @param delimiter Character used as separator.
 * @return Vector of substrings.
 */
std::vector<std::string> split(const std::string& s, char delimiter);

/**
 * @brief Converts a double to a fixed-format string with a given precision.
 *
 * The output uses std::fixed formatting and @p precision digits after
 * the decimal point.
 *
 * @param value     Floating-point value to convert.
 * @param precision Number of digits after the decimal point.
 * @return Formatted string representation of @p value.
 */
std::string doubleToString(double value, int precision);

/**
 * @brief Returns a lowercase copy of the input string.
 *
 * All characters are converted using std::tolower with unsigned char
 * promotion to avoid undefined behavior.
 *
 * @param str Input string.
 * @return Lowercased copy of @p str.
 */
std::string to_lowercase(const std::string& str);

/**
 * @brief Extracts the key set from a std::map into an unordered_set.
 *
 * This utility is useful when one is interested only in the set of keys
 * (e.g. to compare maps by their domains rather than by their values).
 *
 * @tparam T Key type of the map.
 * @tparam U Mapped type of the map.
 * @param map Input std::map.
 * @return An unordered_set containing all keys of @p map.
 */
template<typename T, typename U>
inline std::unordered_set<T> get_keys(const std::map<T, U>& map) {
    std::unordered_set<T> keys;
    for (const std::pair<const  T, U>& item : map) {
        keys.emplace(item.first);
    }
    return keys;
}

/**
 * @brief Extracts the key set from a std::unordered_map into an unordered_set.
 *
 * @tparam T Key type of the map.
 * @tparam U Mapped type of the map.
 * @param map Input std::unordered_map.
 * @return An unordered_set containing all keys of @p map.
 */
template<typename T, typename U>
inline std::unordered_set<T> get_keys(const std::unordered_map<T, U>& map) {
    std::unordered_set<T> keys;
    for (const std::pair<const T, U>& item : map) {
        keys.emplace(item.first);
    }
    return keys;
}

/**
 * @brief Checks whether two std::map instances have the same set of keys.
 *
 * The associated values are ignored; only the key sets are compared.
 *
 * @tparam T Key type of the maps.
 * @tparam U Mapped type of the maps.
 * @param map1 First map.
 * @param map2 Second map.
 * @return true if both maps have identical key sets, false otherwise.
 */
template<typename T, typename U>
inline bool haveSameKeys(const std::map<T, U>& map1, const std::map<T, U>& map2) {
    return get_keys(map1) == get_keys(map2);
}

/**
 * @brief Generic helper to print a value to an output stream.
 *
 * This overload simply forwards the value to the stream using operator<<.
 * It serves as the default implementation for print_value and can be
 * specialized/overloaded for more complex types.
 *
 * @tparam T Value type.
 * @param os    Output stream.
 * @param value Value to print.
 */
template <typename T>
void print_value(std::ostream& os, const T& value) {
    os << value;
}

/**
 * @brief Specialization of print_value for std::shared_ptr.
 *
 * If the pointer is non-null, the pointed-to object is printed. Otherwise
 * the string "nullptr" is written to the stream.
 *
 * @tparam T Pointee type.
 * @param os  Output stream.
 * @param ptr Shared pointer to print.
 */
template <typename T>
void print_value(std::ostream& os, const std::shared_ptr<T> &ptr) {
    if (ptr)
        os << *ptr;
    else
        os << "nullptr";
}

/**
 * @brief Generic stream output operator for map-like containers.
 *
 * This template is enabled only for types Map whose value_type matches
 * std::pair<const key_type, mapped_type>, i.e. true associative maps.
 * The contents are printed in a JSON-like style:
 *
 *   { key1: value1, key2: value2, ... }
 *
 * where values are printed using print_value(), so pointer-like values
 * are handled gracefully.
 *
 * @tparam Map Map-like container type.
 * @param os Output stream.
 * @param m  Map to print.
 * @return Reference to the output stream.
 */
template <typename Map>
std::enable_if_t<
    std::is_same_v<typename Map::value_type,
                   std::pair<const typename Map::key_type, typename Map::mapped_type>>,
    std::ostream&
> operator<<(std::ostream& os, const Map& m) {
    os << "{ ";
    for (const auto& [key, value] : m) {
        os << key << ": ";
        print_value(os, value);
        os << ", ";
    }
    if (!m.empty())
        os.seekp(-2, std::ios_base::end); // remove ", "
    os << " }";
    return os;
}

/**
 * @brief Fills a lookup cache for a function on a finite interval [a, b].
 *
 * The function @p f is evaluated at @p cache_size sample points in the
 * interval [a, b], and the results are stored in @p cache. The first and
 * last points are slightly shifted away from the exact endpoints using
 * nld.epsilon() to avoid potential issues (e.g. singularities) exactly at
 * the boundaries.
 *
 * The intermediate points are placed uniformly between a and b.
 *
 * @tparam T          Value type stored in the cache.
 * @tparam cache_size Size of the cache array.
 * @tparam Func       Callable type; must be invocable as f(double, Args...).
 * @tparam Args       Additional argument types forwarded to @p f.
 *
 * @param f      Function to sample.
 * @param a      Lower bound of the interval.
 * @param b      Upper bound of the interval.
 * @param cache  Preallocated array to fill with sampled values.
 * @param args   Additional arguments forwarded to @p f.
 */
template <typename T, std::size_t cache_size, typename Func, typename... Args>
void fill_cache(Func&& f, double a, double b, std::array<T, cache_size>& cache, Args&&... args) {
    double x_0 = (b - a) * nld.epsilon() + a;
    cache[0] = f(x_0, std::forward<Args>(args)...);
    for (std::size_t i = 1; i < cache_size - 1; ++i) {
        double s = static_cast<double>(i) / static_cast<double>(cache_size - 1);
        double x = (b - a) * s + a;
        cache[i] = f(x, std::forward<Args>(args)...);
    }
    double x_1 = (b - a) * (1. - nld.epsilon()) + a;
    cache[cache_size - 1] = f(x_1, std::forward<Args>(args)...);
}

/**
 * @brief Linearly interpolates a cached function on [a, b].
 *
 * Given a lookup table @p lookup filled with values of some function f at
 * uniformly spaced points in [a, b] (as produced by fill_cache), this
 * routine returns an interpolated value for an input @p x.
 *
 * The input @p x is first clamped to [a, b]. Then:
 *   - it is mapped to a floating-point index in [0, cache_size - 1],
 *   - the two neighboring cache entries are linearly combined according
 *     to the fractional part of the index.
 *
 * @tparam T          Value type stored in the cache and returned.
 * @tparam U          Type of @p x (convertible to double).
 * @tparam cache_size Size of the lookup array.
 *
 * @param x       Input coordinate at which the function is approximated.
 * @param lookup  Cache/lookup table of sampled values.
 * @param a       Lower bound of the sampling interval (default: 0.0).
 * @param b       Upper bound of the sampling interval (default: 1.0).
 * @return Interpolated value at @p x.
 */
template <typename T, typename U, std::size_t cache_size>
T lerp(U x, const std::array<T, cache_size>& lookup, double a=0.0, double b=1.0)
{
    // Clamp s ∈ [0, 1]
    x = std::clamp(static_cast<double>(x), a, b);
    double u = (x - a) / (b - a);
    double scaled = u * (cache_size - 1);
    std::size_t i = static_cast<std::size_t>(std::floor(scaled));
    
    if (i >= cache_size - 1) {
        return lookup[cache_size - 1];
    }

    double t = scaled - static_cast<double>(i);
    return T((1.0 - t)) * lookup[i] + T(t) * lookup[i + 1];
}


#endif // HYPERISO_UTILS_H
