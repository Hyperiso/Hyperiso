#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <optional>
#include <algorithm>
#include <cctype>

// SymbolId<Tag> :  wrapper string-strong-typed
template<class Tag>
class SymbolId {
    std::string name_;
public:
    SymbolId() = default;
    explicit SymbolId(std::string n) : name_(std::move(n)) {}
    const std::string& str() const { return name_; }
};





// ------------------------------------------------------------------
// case-insensitive
// ------------------------------------------------------------------
inline std::string normalize_key(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s) {
        // if (!std::isalnum(c)) continue; // <-- option : ignore space and punctuation
        out.push_back(static_cast<char>(std::tolower(c)));
    }
    return out;
}

template<class Tag>
inline bool operator==(const SymbolId<Tag>& a, const SymbolId<Tag>& b) noexcept {
    return normalize_key(a.str()) == normalize_key(b.str());
}

template<class Tag>
inline bool operator<(const SymbolId<Tag>& a, const SymbolId<Tag>& b) noexcept {
    const auto an = normalize_key(a.str());
    const auto bn = normalize_key(b.str());
    return an < bn;
}

namespace std {
    template<class Tag>
    struct hash<SymbolId<Tag>> {
        size_t operator()(const SymbolId<Tag>& id) const noexcept {
            return std::hash<std::string>{}(normalize_key(id.str()));
        }
    };
}


// without external key
template<class Tag, class ExternalKey, class Hash = std::hash<ExternalKey>>
class DynamicRegistry {
public:
    struct Entry {
        SymbolId<Tag> id;
        std::vector<std::string> aliases;
        std::optional<ExternalKey> ext;
        bool builtin = false;
    };

private:
    std::unordered_map<std::string, Entry> by_name_;
    std::unordered_map<ExternalKey, std::string, Hash> by_ext_;

public:
    bool register_id(const SymbolId<Tag>& id,
                     std::vector<std::string> aliases,
                     std::optional<ExternalKey> ext,
                     bool builtin)
    {
        const std::string canon = normalize_key(id.str());
        Entry e{ id, std::move(aliases), std::move(ext), builtin };

        auto [it, inserted] = by_name_.emplace(canon, e);
        if (!inserted) {
            if (it->second.builtin) return false;   // do not replace builtin
            it->second = e;
        }

        // indexing alias
        for (const auto& a : it->second.aliases) {
            const std::string ak = normalize_key(a);
            auto jt = by_name_.find(ak);
            if (jt == by_name_.end()) {
                by_name_.emplace(ak, it->second);
            } else if (!jt->second.builtin) {
                jt->second = it->second;
            }
        }

        if (it->second.ext) {
            by_ext_[*it->second.ext] = canon;
        }
        return true;
    }

    std::optional<SymbolId<Tag>> find(std::string_view name) const {
        const std::string key = normalize_key(name);
        auto it = by_name_.find(key);
        if (it == by_name_.end()) return std::nullopt;
        return it->second.id;
    }

    std::optional<SymbolId<Tag>> find_by_external(const ExternalKey& k) const {
        auto it = by_ext_.find(k);
        if (it == by_ext_.end()) return std::nullopt;
        auto it2 = by_name_.find(it->second);
        if (it2 == by_name_.end()) return std::nullopt;
        return it2->second.id;
    }

    // key flha from index
    std::optional<ExternalKey> external_of(const SymbolId<Tag>& id) const {
        const std::string key = normalize_key(id.str());
        auto it = by_name_.find(key);
        if (it == by_name_.end() || !it->second.ext) return std::nullopt;
        return it->second.ext;
    }

    bool update_external(const SymbolId<Tag>& id, const ExternalKey& k) {
        const std::string key = normalize_key(id.str());
        auto it = by_name_.find(key);
        if (it == by_name_.end()) return false;
        it->second.ext = k;
        by_ext_[k] = key;
        return true;
    }

    std::vector<SymbolId<Tag>> list_all() const {
        std::vector<SymbolId<Tag>> out;
        std::unordered_map<std::string, bool> seen;
        out.reserve(by_name_.size());
        for (const auto& [_, e] : by_name_) {
            const std::string ck = normalize_key(e.id.str());
            if (!seen.emplace(ck, true).second) continue;
            out.push_back(e.id);
        }
        return out;
    }
};

// spec with external key
template<class Tag>
class DynamicRegistry<Tag, void, void> {
public:
    struct Entry {
        SymbolId<Tag> id;
        std::vector<std::string> aliases;
        bool builtin = false;
    };
private:
    std::unordered_map<std::string, Entry> by_name_;

public:
    bool register_id(const SymbolId<Tag>& id,
                     std::vector<std::string> aliases,
                     bool builtin)
    {
        const std::string canon = normalize_key(id.str());
        Entry e{ id, std::move(aliases), builtin };

        auto [it, inserted] = by_name_.emplace(canon, e);
        if (!inserted) {
            if (it->second.builtin) return false;
            it->second = e;
        }

        for (const auto& a : it->second.aliases) {
            const std::string ak = normalize_key(a);
            auto jt = by_name_.find(ak);
            if (jt == by_name_.end()) {
                by_name_.emplace(ak, it->second);
            } else if (!jt->second.builtin) {
                jt->second = it->second;
            }
        }
        return true;
    }

    std::optional<SymbolId<Tag>> find(std::string_view name) const {
        const std::string key = normalize_key(name);
        auto it = by_name_.find(key);
        if (it == by_name_.end()) return std::nullopt;
        return it->second.id;
    }

    std::vector<SymbolId<Tag>> list_all() const {
        std::vector<SymbolId<Tag>> out;
        std::unordered_map<std::string, bool> seen;
        out.reserve(by_name_.size());
        for (const auto& [_, e] : by_name_) {
            const std::string ck = normalize_key(e.id.str());
            if (!seen.emplace(ck, true).second) continue;
            out.push_back(e.id);
        }
        return out;
    }
};
