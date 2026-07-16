
template<class Tag, class ExternalKey, class Hash>
bool DynamicRegistry<Tag, ExternalKey, Hash>::register_id(const SymbolId<Tag>& id,
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

template<class Tag, class ExternalKey, class Hash>
std::optional<SymbolId<Tag>> DynamicRegistry<Tag, ExternalKey, Hash>::find(std::string_view name) const {
    const std::string key = normalize_key(name);
    auto it = by_name_.find(key);
    if (it == by_name_.end()) return std::nullopt;
    return it->second.id;
}

template<class Tag, class ExternalKey, class Hash>
std::optional<SymbolId<Tag>> DynamicRegistry<Tag, ExternalKey, Hash>::find_by_external(const ExternalKey& k) const {
    auto it = by_ext_.find(k);
    if (it == by_ext_.end()) return std::nullopt;
    auto it2 = by_name_.find(it->second);
    if (it2 == by_name_.end()) return std::nullopt;
    return it2->second.id;
}

template<class Tag, class ExternalKey, class Hash>
std::optional<ExternalKey> DynamicRegistry<Tag, ExternalKey, Hash>::external_of(const SymbolId<Tag>& id) const {
    const std::string key = normalize_key(id.str());
    auto it = by_name_.find(key);
    if (it == by_name_.end() || !it->second.ext) return std::nullopt;
    return it->second.ext;
}

template<class Tag, class ExternalKey, class Hash>
bool DynamicRegistry<Tag, ExternalKey, Hash>::update_external(const SymbolId<Tag>& id, const ExternalKey& k) {
    const std::string key = normalize_key(id.str());
    auto it = by_name_.find(key);
    if (it == by_name_.end()) return false;
    it->second.ext = k;
    by_ext_[k] = key;
    return true;
}

template<class Tag, class ExternalKey, class Hash>
std::vector<SymbolId<Tag>> DynamicRegistry<Tag, ExternalKey, Hash>::list_all() const {
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


template<class Tag>
bool DynamicRegistry<Tag, void, void>::register_id(const SymbolId<Tag>& id,
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

template<class Tag>
std::optional<SymbolId<Tag>> DynamicRegistry<Tag, void, void>::find(std::string_view name) const {
    const std::string key = normalize_key(name);
    auto it = by_name_.find(key);
    if (it == by_name_.end()) return std::nullopt;
    return it->second.id;
}

template<class Tag>
std::vector<SymbolId<Tag>> DynamicRegistry<Tag, void, void>::list_all() const {
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