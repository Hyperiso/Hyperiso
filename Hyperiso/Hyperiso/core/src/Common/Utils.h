#ifndef __HYPERISO_UTILS_H__
#define __HYPERISO_UTILS_H__

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

typedef std::complex<double> complex_t; 

inline std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> parts;
    for (auto &&part : std::views::split(s, delimiter)) {
        parts.emplace_back(std::string_view{&*part.begin(), static_cast<size_t>(std::ranges::distance(part))});
    }
    return parts;
}

inline std::string doubleToString(double value, int precision) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(precision) << value;
	return out.str();
}

inline std::string to_lowercase(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c){ return std::tolower(c); });
    return result;
}

template<typename T, typename U>
inline std::unordered_set<T> get_keys(const std::map<T, U>& map) {
    std::unordered_set<T> keys;
    for (const std::pair<T, U>& item : map) {
        keys.emplace(item.first);
    }
    return keys;
}

template<typename T, typename U>
inline std::unordered_set<T> get_keys(const std::unordered_map<T, U>& map) {
    std::unordered_set<T> keys;
    for (const std::pair<T, U>& item : map) {
        keys.emplace(item.first);
    }
    return keys;
}

template<typename T, typename U>
inline bool haveSameKeys(const std::map<T, U>& map1, const std::map<T, U>& map2) {
    return get_keys(map1) == get_keys(map2);
}

#endif // __HYPERISO_UTILS_H__
