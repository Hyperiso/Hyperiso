#include "Utils.h"

bool ends_with(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> parts;
    std::string::size_type start = 0;

    while (true) {
        auto pos = s.find(delimiter, start);
        if (pos == std::string::npos) {
            parts.emplace_back(s.substr(start));
            break;
        }
        parts.emplace_back(s.substr(start, pos - start));
        start = pos + 1;
    }

    return parts;
}

std::string doubleToString(double value, int precision) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(precision) << value;
	return out.str();
}

std::string to_lowercase(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c){ return std::tolower(c); });
    return result;
}



