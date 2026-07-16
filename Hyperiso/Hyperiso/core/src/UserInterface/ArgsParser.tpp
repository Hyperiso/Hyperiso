#include "ArgsParser.h"

template<typename T>
T ArgParser::convertFromString(const std::string& s) {
    if constexpr (std::is_same_v<T, std::string>) {
        return s;
    } else if constexpr (std::is_same_v<T, int>) {
        size_t pos = 0;
        int v = std::stoi(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("Invalid int: " + s);
        return v;
    } else if constexpr (std::is_same_v<T, double>) {
        size_t pos = 0;
        double v = std::stod(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("Invalid double: " + s);
        return v;
    } else {
        static_assert(!sizeof(T), "Unsupported type for ArgParser::get<T> (use int, double, std::string).");
    }
}

template<typename T>
T ArgParser::get(const std::string& name) const {
    const Argument* argMeta = findArgumentByLongName(name);
    if (!argMeta) {
        throw std::invalid_argument("Unknown argument: " + name);
    }

    if constexpr (std::is_same_v<T, int>) {
        if (argMeta->getType() != ArgType::INT) throw std::invalid_argument("Type mismatch for " + name + " (expected INT)");
    } else if constexpr (std::is_same_v<T, double>) {
        if (argMeta->getType() != ArgType::DOUBLE) throw std::invalid_argument("Type mismatch for " + name + " (expected DOUBLE)");
    } else if constexpr (std::is_same_v<T, std::string>) {
        if (argMeta->getType() != ArgType::STRING) throw std::invalid_argument("Type mismatch for " + name + " (expected STRING)");
    }

    std::string raw = getValue(name);
    return convertFromString<T>(raw);
}

template<typename T>
T ArgParser::getOr(const std::string& name, const T& defaultValue) const {
    if (!exists(name)) return defaultValue;
    return get<T>(name);
}

template<typename T>
std::vector<T> ArgParser::getMany(const std::string& name) const {
    const Argument* argMeta = findArgumentByLongName(name);
    if (!argMeta) {
        throw std::invalid_argument("Unknown argument: " + name);
    }

    if constexpr (std::is_same_v<T, int>) {
        if (argMeta->getType() != ArgType::INT) throw std::invalid_argument("Type mismatch for " + name + " (expected INT)");
    } else if constexpr (std::is_same_v<T, double>) {
        if (argMeta->getType() != ArgType::DOUBLE) throw std::invalid_argument("Type mismatch for " + name + " (expected DOUBLE)");
    } else if constexpr (std::is_same_v<T, std::string>) {
        if (argMeta->getType() != ArgType::STRING) throw std::invalid_argument("Type mismatch for " + name + " (expected STRING)");
    }

    auto raws = getValues(name);
    std::vector<T> out;
    out.reserve(raws.size());
    for (const auto& r : raws) out.push_back(convertFromString<T>(r));
    return out;
}

template<typename T>
std::vector<T> ArgParser::getManyOr(const std::string& name, const std::vector<T>& defaultValues) const {
    if (!exists(name)) return defaultValues;
    return getMany<T>(name);
}
