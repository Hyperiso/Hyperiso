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

// enum class Observables {
//     BR_BS_MUMU,
//     BR_BS_MUMU_UNTAG,
//     BR_BD_MUMU,
//     R_TAU_NU,
//     BR_BU_TAU_NU,
//     ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
//     BR_B_XS_GAMMA,
//     BR_B__D_TAU_NU,
//     A_FB_B__D_TAU_NU,
//     P_TAU_B__D_TAU_NU,
//     R_D,
//     BR_B__DSTAR_TAU_NU,
//     A_FB_B__DSTAR_TAU_NU,
//     P_TAU_B__DSTAR_TAU_NU,
//     P_D_B__DSTAR_TAU_NU,
//     R_DSTAR,
// };

// class ObservableMapper {
// public:
//     static std::string str(Observables obs) {
//         return ObservableMapper::mapping.at(obs);
//     };

//     static Observables enum_elt(std::string name) {
//         return ObservableMapper::inverse_mapping.at(name);
//     };

//     static LhaID flha(Observables obs) {
//         return ObservableMapper::flha_mapping.at(obs);
//     };

//     static Observables enum_elt(LhaID id) {
//         return ObservableMapper::inverse_flha_mapping.at(id);
//     };

//     static std::vector<std::string> get_str() {
//         std::vector<std::string> _;
//         for (auto&& elem : ObservableMapper::mapping) {
//             _.push_back(elem.second);
//         }
//         return _;
//     }

//     static std::vector<Observables> get_enum() {
//         std::vector<Observables> _;
//         for (auto&& elem : ObservableMapper::mapping) {
//             _.push_back(elem.first);
//         }
//         return _;
//     }
// private:
//     static const std::map<Observables, std::string> mapping; 
//     static const std::map<std::string, Observables> inverse_mapping; 
//     static const std::map<Observables, LhaID> flha_mapping; 
//     static const std::map<LhaID, Observables> inverse_flha_mapping; 
// };

// enum class Decays {
//     B__D_l_nu,
//     B__Dstar_l_nu,
//     B__Kstar,
//     B__l_l,
//     B__l_nu,
//     B__Xs,
// };

// class DecayMapper {
// public:
//     static std::vector<Observables> get_observables(Decays decay) {
//         return DecayMapper::obs_mapping.at(decay);
//     }

//     static Decays get_decay(Observables obs) {
//         for (auto &[k, v] : DecayMapper::obs_mapping) {
//             if (std::find(v.begin(), v.end(), obs) != v.end()) {
//                 return k;
//             }
//         }
//         LOG_ERROR("ValueError", "Observable belongs to none of the implemented decays.");
//     }

// private:
//     static const std::map<Decays, std::vector<Observables>> obs_mapping;
// };

// enum class MassType {
//     POLE,
//     MSBAR
// };

// enum class ScaleType {
//     MATCHING,
//     HADRONIC
// };

// class ScaleTypeMapper {
// public:
//     static std::string block(ScaleType type) {
//         return ScaleTypeMapper::block_mapping.at(type);
//     }

// private:
//     static const std::map<ScaleType, std::string> block_mapping; 
// };

// enum class QCDOrder {
//     NONE,
//     LO,
//     NLO,
//     NNLO
// };

// class OrderMapper {

// public:
//     static std::string str(QCDOrder order) {
//         return OrderMapper::mapping.at(order);
//     };

//     static QCDOrder enum_elt(std::string order) {
//         return OrderMapper::inverse_mapping.at(order);
//     };

//     static std::vector<std::string> get_str() {
//         std::vector<std::string> _;
//         for (auto&& elem : OrderMapper::mapping) {
//             _.push_back(elem.second);
//         }
//         return _;
//     }

//     static std::vector<QCDOrder> get_enum() {
//         std::vector<QCDOrder> _;
//         for (auto&& elem : OrderMapper::mapping) {
//             _.push_back(elem.first);
//         }
//         return _;
//     }
// private:
//     static const std::map<QCDOrder, std::string> mapping; 
//     static const std::map<std::string, QCDOrder> inverse_mapping; 
// };

// enum class WCoef {
//     C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2, CBlnu_A, CBlnu_P, C_V1, C_V2, C_S1, C_S2, C_T
// };

// enum class WGroup {
//     B, 
//     BPrime, 
//     BScalar,
//     Blnu,
//     BCLNU,
// };

// enum class BWilsonBasis {
//     STANDARD, 
//     TRADITIONAL
// };

// enum class ContributionType {
//     SM, 
//     BSM,
//     TOTAL
// };


// class GroupMapper {
// public:
//     static std::string str(WGroup group) {
//         return GroupMapper::mapping.at(group);
//     };

//     static WGroup enum_elt(std::string group) {
//         return GroupMapper::inverse_mapping.at(group);
//     };

//     static std::vector<std::string> get_str() {
//         std::vector<std::string> _;
//         for (auto&& elem : GroupMapper::mapping) {
//             _.push_back(elem.second);
//         }
//         return _;
//     }

//     static std::vector<WGroup> get_enum() {
//         std::vector<WGroup> _;
//         for (auto&& elem : GroupMapper::mapping) {
//             _.push_back(elem.first);
//         }
//         return _;
//     }
// private:
//     static const std::map<WGroup, std::string> mapping; 
//     static const std::map<std::string, WGroup> inverse_mapping; 
// };

// class WCoefMapper {

// public:
//     static std::string str(WCoef coef) {
//         return WCoefMapper::mapping.at(coef);
//     };

//     static WCoef enum_elt(std::string coef) {
//         return WCoefMapper::inverse_mapping.at(coef);
//     };

//     static std::pair<int, int> flha_base(WCoef coef) {
//         return WCoefMapper::flha_mapping.at(coef);
//     };

//     static LhaID flha_full(WCoef coef, QCDOrder order, ContributionType type) {
//         // LOG_INFO("Attention");
//         auto base_id = WCoefMapper::flha_mapping.at(coef);
//         // LOG_INFO("Paf");
//         return LhaID{base_id.first, base_id.second, static_cast<int>(order) - 1, static_cast<int>(type)};
//     };

//     static WCoef from_flha(int content, int structure) {
//         if (WCoefMapper::inverse_flha_mapping.contains({content, structure})) {
//             return WCoefMapper::inverse_flha_mapping.at({content, structure});
//         }
//         LOG_ERROR("General", "Wilson coefficient with ID", content, structure, "is not supported");
//     };

//     static std::vector<WCoef> get_group(WGroup group) {
//         switch (group) {
//             case WGroup::B:
//                 return B_group;
//             case WGroup::BPrime:
//                 return B_prime_group;
//             case WGroup::BScalar:
//                 return B_scalar_group;
//             case WGroup::Blnu:
//                 return B_lnu_group;
//             case WGroup::BCLNU:
//                 return b_clnu_group;
//             default:
//                 LOG_ERROR("Invalid WGroup", "get_group function couldn't find your group");
//         }
//     }

//     static size_t n_wilsons() {
//         return WCoefMapper::mapping.size();
//     }

//     static std::vector<std::string> get_str() {
//         std::vector<std::string> _;
//         for (auto&& elem : WCoefMapper::mapping) {
//             _.push_back(elem.second);
//         }
//         return _;
//     }

//     static std::vector<WCoef> get_enum() {
//         std::vector<WCoef> _;
//         for (auto&& elem : WCoefMapper::mapping) {
//             _.push_back(elem.first);
//         }
//         return _;
//     }
// private:
//     static const std::vector<WCoef> B_group;
//     static const std::vector<WCoef> B_prime_group;
//     static const std::vector<WCoef> B_scalar_group;
//     static const std::vector<WCoef> B_lnu_group;
//     static const std::vector<WCoef> b_clnu_group;
//     static const std::map<WCoef, std::string> mapping; 
//     static const std::map<std::string, WCoef> inverse_mapping; 
//     static const std::map<WCoef, std::pair<int, int>> flha_mapping; 
//     static const std::map<std::pair<int, int>, WCoef> inverse_flha_mapping; 
// };

/* !!!! Do not change the order of the first 4 entries !!!! */

// enum class ParameterType {
//     SM,
//     BSM,
//     FLAVOR,
//     WILSON,
//     DECAY,
//     PASSTHROUGH,
//     OBSERVABLE,
// };

// class ParameterTypeMapper {
// public:
//     static std::string str(ParameterType type) {
//         return ParameterTypeMapper::mapping.at(type);
//     };

//     static ParameterType enum_elt(std::string type) {
//         return ParameterTypeMapper::inverse_mapping.at(type);
//     };

//     static std::vector<std::string> get_str() {
//         std::vector<std::string> _;
//         for (auto&& elem : ParameterTypeMapper::mapping) {
//             _.push_back(elem.second);
//         }
//         return _;
//     }

//     static std::vector<ParameterType> get_enum() {
//         std::vector<ParameterType> _;
//         for (auto&& elem : ParameterTypeMapper::mapping) {
//             _.push_back(elem.first);
//         }
//         return _;
//     }

// private:
//     static const std::map<ParameterType, std::string> mapping; 
//     static const std::map<std::string, ParameterType> inverse_mapping; 
// };

// enum class Model {
//     SM,
//     SUSY,
//     THDM,
//     CUSTOM
// };

// class ModelMapper {
// public:
//     static std::string str(Model model) {
//         return ModelMapper::mapping.at(model);
//     };

//     static Model enum_elt(std::string model) {
//         return ModelMapper::inverse_mapping.at(model);
//     };

//     static std::vector<std::string> get_str() {
//         std::vector<std::string> _;
//         for (auto&& elem : ModelMapper::mapping) {
//             _.push_back(elem.second);
//         }
//         return _;
//     }

//     static std::vector<Model> get_enum() {
//         std::vector<Model> _;
//         for (auto&& elem : ModelMapper::mapping) {
//             _.push_back(elem.first);
//         }
//         return _;
//     }

// private:
//     static const std::map<Model, std::string> mapping; 
//     static const std::map<std::string, Model> inverse_mapping; 
// };

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
            return std::hash<BlockName>{}(p.block) ^ std::hash<LhaID>{}(p.code);
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const ParamId& pid) {
    os << pid.block << ":" << pid.code;
    return os;
};


class DependenciesHelper {
public:
    static std::unordered_set<ParamId> get_allowed_parameters(Observables id);
    static bool is_param_allowed(Observables id, ParamId pid);

private:
    static const std::map<Observables, std::unordered_set<ParamId>> dep_lists;
};

class LhaParamsHelper {
public:
    static std::vector<std::vector<long>> get_minimal_content(const BlockName& block_name);

private:
    static const std::map<BlockName, std::vector<std::vector<long>>> minimal_blocks;
};

#endif // __GENERAL_H__