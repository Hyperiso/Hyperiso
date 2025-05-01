#ifndef ENUM_MAPPER_H
#define ENUM_MAPPER_H
#include <map>
#include <string>
#include <vector>

template <typename Enum, typename Derived>
class EnumMapperBase {
public:
    static std::string str(Enum e) {
        return Derived::mapping().at(e);
    }

    static Enum enum_elt(const std::string& name) {
        return Derived::inverse_mapping().at(name);
    }

    static std::vector<std::string> get_str() {
        std::vector<std::string> out;
        for (const auto& [key, val] : Derived::mapping()) {
            out.push_back(val);
        }
        return out;
    }

    static std::vector<Enum> get_enum() {
        std::vector<Enum> out;
        for (const auto& [key, val] : Derived::mapping()) {
            out.push_back(key);
        }
        return out;
    }
};

template <typename K, typename V>
std::map<V, K> invert_map(const std::map<K, V>& original) {
    std::map<V, K> inv;
    for (const auto& [k, v] : original) {
        inv[v] = k;
    }
    return inv;
}

#endif