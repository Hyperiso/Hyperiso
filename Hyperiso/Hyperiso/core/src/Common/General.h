#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <string>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <initializer_list>
#include <map>
#include <vector>
#include <set>
#include <ranges>
#include <concepts>
#include <variant>
#include "Logger.h"
#include "Utils.h"
#include "GeneralEnum.h"

namespace fs = std::filesystem;

/**
 * @struct LhaID
 * @brief Represents an identifier of a LHA element, possibly containing several sub-ids 
 */
struct LhaID {
    std::vector<long> parts;     /**< Collection of sub-ids. */

    /**
     * @brief Constructs a LhaID with specified sub-ids
     * @param parts List of sub-ids of the element
     */
    template<typename... Args>
    requires (std::convertible_to<Args, long> && ...)
    LhaID(Args... sub_ids) : parts({static_cast<long>(sub_ids)...}) {}

    /**
     * @brief Constructs a LhaID from a string of _-separated integer values
     * @param parts String of sub-ids of the element separated by an _
     */
    LhaID(const std::string& parts);

    /**
     * @brief Constructs a LhaID with a single identifier
     * @param sub_ids Sub-identifiers of the element
     */
    LhaID(const std::vector<long>& sub_ids) : parts(std::move(sub_ids)) {}

    LhaID(std::initializer_list<long> sub_ids) : parts(sub_ids) {}

    /**
     * @brief Constructs a LhaID with a single identifier
     * @param id Identifier of the element
     */
    LhaID(long id) : parts({id}) {}

    std::string to_string() const;

    std::vector<long> get_parts() const { return parts; }
    
    /**
     * @brief Allows for implicit conversion of a trivial LhaID to an integer 
     */
    operator long() const {
        if (this->parts.size() > 1) {
            LOG_WARN("Casting nontrivial LhaID to int discards information.");
            for (auto& part : this->parts){
                std::cout << part<< std::endl;
            }
        }
        
        return this->parts.at(0);
    };

    inline friend bool operator==(const LhaID& lhs, const LhaID& rhs) { return lhs.parts == rhs.parts; };
    inline friend bool operator!=(const LhaID& lhs, const LhaID& rhs) { return !(lhs == rhs); };
    inline friend bool operator<(const LhaID& lhs, const LhaID& rhs) { return (lhs.parts <=> rhs.parts) == std::weak_ordering::less; };

    friend std::ostream& operator<<(std::ostream&, const LhaID&);
};

namespace std {
    template <>
    struct hash<LhaID> {
        std::size_t operator()(const LhaID& p) const noexcept {
            return std::hash<std::string>{}(p.to_string());
        }
    };
}

class BlockName {
private:
    std::unordered_set<std::string> block_names;

public:
    BlockName() = default;

    BlockName(const std::string& name) {
        block_names.insert(name);
    }

    BlockName(const char* name) {
        block_names.insert(std::string(name));
    }
    
    BlockName(std::initializer_list<std::string> names) : block_names(names) {}

    BlockName(const std::unordered_set<std::string>& names) : block_names(names) {}

    std::unordered_set<std::string> get_alias() const{
        return block_names;
    }
    operator std::string() const {
        if (block_names.size() > 1) {
            LOG_WARN("Casting BlockName with multiple aliases to string discards information.");
            for (const auto& name : block_names) {
                std::cerr << name << std::endl;
            }
        }
        return block_names.empty() ? "" : *block_names.begin();
    }

    bool hasAlias(const std::string& alias) const {
        return block_names.find(alias) != block_names.end();
    }

    std::string to_string() const {
        return operator std::string();
    }

    bool operator==(const BlockName& other) const {
        for (const auto& name : block_names) {
            if (other.hasAlias(name)) return true;
        }
        return false;
    }

    bool operator!=(const BlockName& other) const {
        return !(*this == other);
    }

    bool operator==(const std::string& name) const {
        return hasAlias(name);
    }

    bool operator==(const char* name) const {
        return hasAlias(std::string(name));
    }

    bool operator!=(const std::string& name) const {
        return !hasAlias(name);
    }

    BlockName& addAlias(const std::string& alias) {
        block_names.insert(alias);
        return *this;
    }

    friend BlockName operator+(const std::string& lhs, const BlockName& rhs) {
        std::unordered_set<std::string> combined;
        for (const auto& name : rhs.block_names) {
            combined.insert(lhs + name);
        }
        return BlockName{combined};
    }

    friend std::ostream& operator<<(std::ostream& os, const BlockName& block) {
        bool first = true;
        for (const auto& name : block.block_names) {
            if (!first) os << "/";
            os << name;
            first = false;
        }
        return os;
    }

    void to_upper() {
        std::unordered_set<std::string> upper_names;
        for (auto name : block_names) {
            std::transform(name.begin(), name.end(), name.begin(), ::toupper);
            upper_names.insert(name);
        }
        block_names = std::move(upper_names);
    }

    bool operator<(const BlockName& other) const {
        std::set<std::string> lhs_sorted(block_names.begin(), block_names.end());
        std::set<std::string> rhs_sorted(other.block_names.begin(), other.block_names.end());
        return lhs_sorted < rhs_sorted;
    }
};

inline bool operator==(const std::string& lhs, const BlockName& rhs) {
    return rhs == lhs;
}

inline bool operator!=(const std::string& lhs, const BlockName& rhs) {
    return !(rhs == lhs);
}


namespace std {
    template <>
    struct hash<BlockName> {
        std::size_t operator()(const BlockName& p) const noexcept {
            size_t h = 0;
            for (const auto& name : p.get_alias()) {
                h ^= std::hash<std::string>{}(name) + 0x9e3779b9 + (h << 6) + (h >> 2); // boost-style hash combine
            }
            return h;
        }
    };
}


struct ParamId {
    std::optional<ParameterType> type;
    BlockName block;
    LhaID code;

    ParamId() : block("NULL"), code(0) {}
    ParamId(const BlockName& block, const LhaID& code) : block(block), code(code) {}
    ParamId(ParameterType type, const BlockName& block, const LhaID& code) : type(type), block(block), code(code) {}

    void set_parameter_type(ParameterType type) { this->type = type; }

    inline friend bool operator==(const ParamId& lhs, const ParamId& rhs) { 
        return lhs.type == rhs.type && lhs.block == rhs.block && lhs.code == rhs.code;
    };

    bool operator<(const ParamId& other) const {
        if (type != other.type) return type < other.type;
        if (block != other.block) return block < other.block;
        return code < other.code;
    }
};

namespace std {
    template <>
    struct hash<ParamId> {
        std::size_t operator()(const ParamId& p) const noexcept {
            std::size_t h1 = std::hash<BlockName>{}(p.block);
            std::size_t h2 = std::hash<LhaID>{}(p.code);
            std::size_t h3 = p.type ? std::hash<int>{}(static_cast<int>(*p.type)) : 0;
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const ParamId& pid) {
    os << pid.block << ":" << pid.code;
    return os;
};




class LhaParamsHelper {
public:
    static std::vector<std::vector<long>> get_minimal_content(const BlockName& block_name);

private:
    static const std::map<BlockName, std::vector<std::vector<long>>> minimal_blocks;
};

#endif // __GENERAL_H__