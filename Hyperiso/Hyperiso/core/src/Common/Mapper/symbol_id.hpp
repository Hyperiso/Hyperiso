#ifndef SYMBOL_ID_HPP
#define SYMBOL_ID_HPP

#include <string>
#include <string_view>
#include <algorithm>
/**
 * @file SymbolId.h
 * @brief Lightweight strongly typed wrapper around string identifiers.
 * !! NOT use for the moment !!
 */

/**
 * @brief Returns an uppercase copy of the input string.
 */
inline std::string to_upper_copy(std::string_view s) {
    std::string out(s.begin(), s.end());
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return out;
}

/**
 * @file SymbolId.h
 * @brief Lightweight strongly typed wrapper around string identifiers.
 */

/**
 * @brief Returns an uppercase copy of the input string.
 */
template <class Tag>
class SymbolId {
public:
    /// Constructs an empty identifier.
    SymbolId() = default;

    /// Constructs an identifier from a string.
    explicit SymbolId(std::string name) : name_(std::move(name)) {}

    /// Returns the underlying string.
    const std::string& str() const { return name_; }

    /// Explicit conversion to const std::string&.
    explicit operator const std::string&() const { return name_; }

    /// Returns true if the identifier is empty.
    bool empty() const { return name_.empty(); }

    friend bool operator==(const SymbolId& a, const SymbolId& b) { return a.name_ == b.name_; }
    friend bool operator!=(const SymbolId& a, const SymbolId& b) { return !(a==b); }
    friend bool operator<(const SymbolId& a, const SymbolId& b)  { return a.name_ < b.name_; }

private:
    std::string name_;
};

/**
 * @brief std::hash specialization for SymbolId<Tag>.
 */
namespace std {
    template<class Tag>
    struct hash<SymbolId<Tag>> {
        std::size_t operator()(const SymbolId<Tag>& id) const noexcept {
            return std::hash<std::string>{}(id.str());
        }
    };
}

#endif