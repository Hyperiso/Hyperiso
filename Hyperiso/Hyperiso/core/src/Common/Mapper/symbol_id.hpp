#pragma once
#include <string>
#include <string_view>
#include <algorithm>

inline std::string to_upper_copy(std::string_view s) {
    std::string out(s.begin(), s.end());
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return out;
}

template <class Tag>
class SymbolId {
public:
    SymbolId() = default;
    explicit SymbolId(std::string name) : name_(std::move(name)) {}

    const std::string& str() const { return name_; }
    explicit operator const std::string&() const { return name_; }
    bool empty() const { return name_.empty(); }

    friend bool operator==(const SymbolId& a, const SymbolId& b) { return a.name_ == b.name_; }
    friend bool operator!=(const SymbolId& a, const SymbolId& b) { return !(a==b); }
    friend bool operator<(const SymbolId& a, const SymbolId& b)  { return a.name_ < b.name_; }

private:
    std::string name_;
};

namespace std {
    template<class Tag>
    struct hash<SymbolId<Tag>> {
        std::size_t operator()(const SymbolId<Tag>& id) const noexcept {
            return std::hash<std::string>{}(id.str());
        }
    };
}
