#include "BlockName.h"

BlockName::BlockName(const std::string& name) {
    block_names.insert(name);
}

BlockName::BlockName(const char* name) {
    block_names.insert(std::string(name));
}

BlockName::BlockName(std::initializer_list<std::string> names)
: block_names(names.begin(), names.end()) {}

BlockName::BlockName(const std::unordered_set<std::string>& names)
: block_names(names) {}

std::unordered_set<std::string> BlockName::get_alias() const {
    return block_names;
}

BlockName::operator std::string() const {
    if (block_names.size() > 1) {
        LOG_TRACE("Casting BlockName with multiple aliases to string discards information.");
        for (const auto& name : block_names) LOG_TRACE(name);
    }
    return block_names.empty() ? "" : *block_names.begin();
}

bool BlockName::hasAlias(const std::string& alias) const {
    return block_names.find(alias) != block_names.end();
}

std::string BlockName::to_string() const {
    return operator std::string();
}

bool BlockName::operator==(const BlockName& other) const {
    // intersection non vide
    for (const auto& name : block_names) {
        if (other.hasAlias(name)) return true;
    }
    return false;
}

bool BlockName::operator!=(const BlockName& other) const {
    return !(*this == other);
}

bool BlockName::operator==(const std::string& name) const {
    return hasAlias(name);
}

bool BlockName::operator==(const char* name) const {
    return hasAlias(std::string(name));
}

bool BlockName::operator!=(const std::string& name) const {
    return !hasAlias(name);
}

BlockName& BlockName::addAlias(const std::string& alias) {
    block_names.insert(alias);
    return *this;
}

BlockName operator+(const std::string& lhs, const BlockName& rhs) {
    std::unordered_set<std::string> combined;
    for (const auto& name : rhs.block_names) {
        combined.insert(lhs + name);
    }
    return BlockName{combined};
}

std::ostream& operator<<(std::ostream& os, const BlockName& block) {
    bool first = true;
    for (const auto& name : block.block_names) {
        if (!first) os << "/";
        os << name;
        first = false;
    }
    return os;
}

void BlockName::to_upper() {
    std::unordered_set<std::string> upper_names;
    upper_names.reserve(block_names.size());
    for (auto name : block_names) {
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        upper_names.insert(std::move(name));
    }
    block_names = std::move(upper_names);
}

bool BlockName::operator<(const BlockName& other) const {
    // ordre stable : compare lexicographiquement les sets triés
    std::set<std::string> lhs_sorted(block_names.begin(), block_names.end());
    std::set<std::string> rhs_sorted(other.block_names.begin(), other.block_names.end());
    return lhs_sorted < rhs_sorted;
}
