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
#include <memory>

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

template <typename T>
void print_value(std::ostream& os, const T& value) {
    os << value;
}

template <typename T>
void print_value(std::ostream& os, const std::shared_ptr<T> &ptr) {
    if (ptr)
        os << *ptr;
    else
        os << "nullptr";
}

template <typename Map>
std::enable_if_t<
    std::is_same_v<typename Map::value_type,
                   std::pair<const typename Map::key_type, typename Map::mapped_type>>,
    std::ostream&
>
operator<<(std::ostream& os, const Map& m) {
    os << "{ ";
    for (const auto& [key, value] : m) {
        os << key << ": ";
        print_value(os, value);
        os << ", ";
    }
    if (!m.empty())
        os.seekp(-2, std::ios_base::end); // retire ", "
    os << " }";
    return os;
}




#endif // __HYPERISO_UTILS_H__
