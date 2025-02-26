#ifndef __HYPERISO_UTILS_H__
#define __HYPERISO_UTILS_H__

#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

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

#endif // __HYPERISO_UTILS_H__
