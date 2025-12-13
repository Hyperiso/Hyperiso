// #include "BlockName.h"


// // BlockName::BlockName(const std::string& name) {
// //     canonical_name = name;
// //     block_names.insert(name);
// // }

// // BlockName::BlockName(const char* name) : BlockName(std::string(name)) {}

// // BlockName::BlockName(const std::string& name)
// // : block_names(split_aliases(name)) {}

// // BlockName::BlockName(const char* name)
// // : BlockName(std::string(name)) {}

// // BlockName::BlockName(const std::string& name) {
// //     block_names.insert(name);
// // }

// // BlockName::BlockName(const std::string& name) {
// //     auto up = upper_copy(name);
// //     // si "A/B/C" on split, primary = premier token (stable)
// //     auto sp = split_aliases(up);
// //     if (!sp.empty()) {
// //         // primary = premier token dans la string originale (pas dans le set)
// //         // on le re-extrait proprement :
// //         auto pos = up.find('/');
// //         primary = (pos == std::string::npos) ? up : up.substr(0, pos);
// //         block_names = std::move(sp);
// //     } else {
// //         primary = up;
// //         block_names.insert(up);
// //     }
// // }

// BlockName::BlockName(const std::string& name) {
//     // split "YU/UCOUPL/YL" -> {"YU","UCOUPL","YL"}
//     auto aliases = split_slash_aliases(name);
//     block_names.insert(aliases.begin(), aliases.end());
// }

// // BlockName::BlockName(const char* name) {
// //     block_names.insert(std::string(name));
// // }

// // BlockName::BlockName(const char* name) : BlockName(std::string(name)) {}

// BlockName::BlockName(const char* name) : BlockName(std::string(name ? name : "")) {}

// // BlockName::BlockName(std::initializer_list<std::string> names) : block_names(names) {
// //     if (!names.size()) return;
// //     canonical_name = *names.begin();
// // }

// // BlockName::BlockName(std::initializer_list<std::string> names) {
// //     if (names.size() == 0) return;
// //     // primary = premier de la liste (stable)
// //     primary = upper_copy(*names.begin());
// //     for (auto& n : names) block_names.insert(upper_copy(n));
// //     block_names.insert(primary);
// // }


// // // BlockName::BlockName(const std::unordered_set<std::string>& names) : block_names(names) {
// // //     canonical_name = block_names.empty() ? "" : *block_names.begin();
// // // }

// // BlockName::BlockName(const std::unordered_set<std::string>& names) {
// //     if (names.empty()) return;
// //     // pas d'ordre dans unordered_set => on choisit un primary stable :
// //     // on prend le plus petit LEXI, mais IMPORTANT : ça ne sert qu'à stabiliser
// //     // un BlockName "construit depuis un set". En pratique toi tu construis
// //     // presque toujours depuis string / init_list.
// //     auto it = std::min_element(names.begin(), names.end());
// //     primary = upper_copy(*it);
// //     for (auto& n : names) block_names.insert(upper_copy(n));
// //     block_names.insert(primary);
// // }

// BlockName::BlockName(std::initializer_list<std::string> names) {
//     for (const auto& s : names) {
//         auto aliases = split_slash_aliases(s);
//         block_names.insert(aliases.begin(), aliases.end());
//     }
// }

// BlockName::BlockName(const std::unordered_set<std::string>& names) {
//     for (const auto& s : names) {
//         auto aliases = split_slash_aliases(s);
//         block_names.insert(aliases.begin(), aliases.end());
//     }
// }

// const std::string& BlockName::canonical() const { return primary; }

// std::unordered_set<std::string> BlockName::get_alias() const{
//     return block_names;
// }

// // BlockName::operator std::string() const {
// //     if (block_names.size() > 1) {
// //         LOG_WARN("Casting BlockName with multiple aliases to string discards information.");
// //         for (const auto& name : block_names) {
// //             std::cerr << name << std::endl;
// //         }
// //     }
// //     return block_names.empty() ? "" : *block_names.begin();
// // }

// // BlockName::operator std::string() const { return primary; }

// BlockName::operator std::string() const {
//     if (block_names.size() > 1) {
//         LOG_WARN("Casting BlockName with multiple aliases to string discards information.");
//     }
//     return block_names.empty() ? "" : *block_names.begin();
// }

// // bool BlockName::hasAlias(const std::string& alias) const {
// //     return block_names.find(alias) != block_names.end();
// // }

// // bool BlockName::hasAlias(const std::string& a) const {
// //     return block_names.contains(upper_copy(a));
// // }

// bool BlockName::hasAlias(const std::string& alias) const {
//     return block_names.find(alias) != block_names.end();
// }

// // std::string BlockName::to_string() const { return primary; }

// std::string BlockName::to_string() const {
//     return operator std::string();
// }
// // std::string BlockName::to_string() const {
// //     return operator std::string();
// // }

// // bool BlockName::operator==(const BlockName& other) const {
// //     for (const auto& name : block_names) {
// //         if (other.hasAlias(name)) return true;
// //     }
// //     return false;
// // }

// // bool BlockName::operator==(const BlockName& other) const { return primary == other.primary; }

// bool BlockName::operator==(const BlockName& other) const {
//     // itère sur le plus petit set pour être un peu plus efficace
//     if (block_names.size() <= other.block_names.size()) {
//         for (const auto& a : block_names) {
//             if (other.block_names.find(a) != other.block_names.end()) return true;
//         }
//     } else {
//         for (const auto& a : other.block_names) {
//             if (block_names.find(a) != block_names.end()) return true;
//         }
//     }
//     return false;
// }

// // bool BlockName::operator==(const BlockName& other) const {
// //     return canonical() == other.canonical();
// // }

// // bool BlockName::operator!=(const BlockName& other) const {
// //     return !(*this == other);
// // }

// // bool  BlockName::operator!=(const BlockName& other) const { return !(*this == other); }

// bool BlockName::operator!=(const BlockName& other) const {
//     return !(*this == other);
// }

// // bool BlockName::operator==(const std::string& name) const {
// //     return hasAlias(name);
// // }

// // bool BlockName::operator==(const char* name) const {
// //     return hasAlias(name);
// // }

// // bool BlockName::operator!=(const std::string& name) const {
// //     return !hasAlias(name);
// // }

// bool BlockName::operator==(const std::string& name) const {
//     // note: ici on NE split PAS le string, on teste l'alias exact
//     // (tes tests attendent b != "mass" après to_upper)
//     return hasAlias(name);
// }

// bool BlockName::operator==(const char* name) const {
//     return hasAlias(std::string(name ? name : ""));
// }

// bool BlockName::operator!=(const std::string& name) const {
//     return !hasAlias(name);
// }

// // BlockName& BlockName::addAlias(const std::string& alias) {
// //     block_names.insert(alias);
// //     return *this;
// // }

// // BlockName& BlockName::addAlias(const std::string& a) {
// //     block_names.insert(upper_copy(a));
// //     return *this;
// // }

// BlockName& BlockName::addAlias(const std::string& alias) {
//     auto aliases = split_slash_aliases(alias);
//     block_names.insert(aliases.begin(), aliases.end());
//     return *this;
// }

// BlockName& BlockName::addAliases(const std::unordered_set<std::string>& as) {
//     for (auto& a : as) block_names.insert(upper_copy(a));
//     return *this;
// }

// // BlockName operator+(const std::string& lhs, const BlockName& rhs) {
// //     std::unordered_set<std::string> combined;
// //     for (const auto& name : rhs.block_names) combined.insert(lhs + name);

// //     BlockName out{combined};
// //     out.canonical_name = lhs + rhs.canonical();
// //     return out;
// // }

// // BlockName operator+(const std::string& lhs, const BlockName& rhs) {
// //     std::unordered_set<std::string> combined;
// //     for (auto& a : rhs.block_names) combined.insert(lhs + a);
// //     return BlockName(combined);
// // }

// BlockName operator+(const std::string& lhs, const BlockName& rhs) {
//     std::unordered_set<std::string> combined;
//     for (const auto& name : rhs.block_names) {
//         combined.insert(lhs + name);
//     }
//     return BlockName{combined};
// }

// std::ostream& operator<<(std::ostream& os, const BlockName& block) {
//     bool first = true;
//     for (const auto& name : block.block_names) {
//         if (!first) os << "/";
//         os << name;
//         first = false;
//     }
//     return os;
// }

// // void BlockName::to_upper() {
// //     std::unordered_set<std::string> upper_names;
// //     for (auto name : block_names) {
// //         std::transform(name.begin(), name.end(), name.begin(), ::toupper);
// //         upper_names.insert(name);
// //     }
// //     block_names = std::move(upper_names);
// // }

// // void BlockName::to_upper() {
// //     primary = upper_copy(primary);
// //     std::unordered_set<std::string> up;
// //     for (auto& a : block_names) up.insert(upper_copy(a));
// //     block_names = std::move(up);
// //     block_names.insert(primary);
// // }

// void BlockName::to_upper() {
//     std::unordered_set<std::string> upper_names;
//     upper_names.reserve(block_names.size());
//     for (auto name : block_names) {
//         std::transform(name.begin(), name.end(), name.begin(), ::toupper);
//         upper_names.insert(std::move(name));
//     }
//     block_names = std::move(upper_names);
// }

// // bool BlockName::operator<(const BlockName& other) const {
// //     std::set<std::string> lhs_sorted(block_names.begin(), block_names.end());
// //     std::set<std::string> rhs_sorted(other.block_names.begin(), other.block_names.end());
// //     return lhs_sorted < rhs_sorted;
// // }

// // bool BlockName::operator<(const BlockName& other) const {
// //     // pour std::map dans LhaParser :contentReference[oaicite:2]{index=2}
// //     return primary < other.primary;
// // }

// bool BlockName::operator<(const BlockName& other) const {
//     std::set<std::string> lhs_sorted(block_names.begin(), block_names.end());
//     std::set<std::string> rhs_sorted(other.block_names.begin(), other.block_names.end());
//     return lhs_sorted < rhs_sorted;
// }
















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
        LOG_WARN("Casting BlockName with multiple aliases to string discards information.");
        for (const auto& name : block_names) std::cerr << name << std::endl;
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
