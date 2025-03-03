#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <string>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <map>
#include <vector>
#include <ranges>
#include <concepts>
#include <variant>
#include "Logger.h"
#include "Utils.h"

enum class Observables {
    BR_BS_MUMU,
    BR_BS_MUMU_UNTAG,
    BR_BD_MUMU,
    R_TAU_NU,
    BR_BU_TAU_NU,
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
    BR_B_XS_GAMMA,
    BR_B__D_TAU_NU,
    A_FB_B__D_TAU_NU,
    P_TAU_B__D_TAU_NU,
    R_D,BR_B__DSTAR_TAU_NU,
    A_FB_B__DSTAR_TAU_NU,
    P_TAU_B__DSTAR_TAU_NU,
    P_D_B__DSTAR_TAU_NU,
    R_DSTAR,
};

class ObservableMapper {
public:
    static std::string str(Observables obs) {
        return ObservableMapper::mapping.at(obs);
    };

    static Observables enum_elt(std::string name) {
        return ObservableMapper::inverse_mapping.at(name);
    };

    static std::vector<std::string> get_str() {
        std::vector<std::string> _;
        for (auto&& elem : ObservableMapper::mapping) {
            _.push_back(elem.second);
        }
        return _;
    }

    static std::vector<Observables> get_enum() {
        std::vector<Observables> _;
        for (auto&& elem : ObservableMapper::mapping) {
            _.push_back(elem.first);
        }
        return _;
    }
private:
    static const std::map<Observables, std::string> mapping; 
    static const std::map<std::string, Observables> inverse_mapping; 
};

enum class Decays {
    B__D_l_nu,
    B__Dstar_l_nu,
    B__Kstar,
    B__l_l,
    B__l_nu,
    B__Xs,
};

class DecayMapper {
public:
    static std::vector<Observables> get_observables(Decays decay) {
        return DecayMapper::obs_mapping.at(decay);
    }

    static Decays get_decay(Observables obs) {
        for (auto &[k, v] : DecayMapper::obs_mapping) {
            if (std::find(v.begin(), v.end(), obs) != v.end()) {
                return k;
            }
        }
        LOG_ERROR("ValueError", "Observable belongs to none of the implemented decays.");
    }

private:
    static const std::map<Decays, std::vector<Observables>> obs_mapping;
};

enum class QCDOrder {
    NONE,
    LO,
    NLO,
    NNLO,
};

class OrderMapper {

public:
    static std::string str(QCDOrder order) {
        return OrderMapper::mapping.at(order);
    };

    static QCDOrder enum_elt(std::string order) {
        return OrderMapper::inverse_mapping.at(order);
    };

    static std::vector<std::string> get_str() {
        std::vector<std::string> _;
        for (auto&& elem : OrderMapper::mapping) {
            _.push_back(elem.second);
        }
        return _;
    }

    static std::vector<QCDOrder> get_enum() {
        std::vector<QCDOrder> _;
        for (auto&& elem : OrderMapper::mapping) {
            _.push_back(elem.first);
        }
        return _;
    }
private:
    static const std::map<QCDOrder, std::string> mapping; 
    static const std::map<std::string, QCDOrder> inverse_mapping; 
};

enum class WCoef {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2, CBlnu_A, CBlnu_P, C_V1, C_V2, C_S1, C_S2, C_T
};

enum class WGroup {
    B, 
    BPrime, 
    BScalar,
    Blnu,
    BCLNU,
};

enum class BWilsonBasis {
    STANDARD, 
    TRADITIONAL
};

class GroupMapper {
public:
    static std::string str(WGroup group) {
        return GroupMapper::mapping.at(group);
    };

    static WGroup enum_elt(std::string group) {
        return GroupMapper::inverse_mapping.at(group);
    };

    static std::vector<std::string> get_str() {
        std::vector<std::string> _;
        for (auto&& elem : GroupMapper::mapping) {
            _.push_back(elem.second);
        }
        return _;
    }

    static std::vector<WGroup> get_enum() {
        std::vector<WGroup> _;
        for (auto&& elem : GroupMapper::mapping) {
            _.push_back(elem.first);
        }
        return _;
    }
private:
    static const std::map<WGroup, std::string> mapping; 
    static const std::map<std::string, WGroup> inverse_mapping; 
};

class WCoefMapper {

public:
    static std::string str(WCoef coef) {
        return WCoefMapper::mapping.at(coef);
    };

    static WCoef enum_elt(std::string coef) {
        return WCoefMapper::inverse_mapping.at(coef);
    };

    static std::pair<int, int> flha(WCoef coef) {
        return WCoefMapper::flha_mapping.at(coef);
    };

    static WCoef from_flha(int content, int structure) {
        if (WCoefMapper::inverse_flha_mapping.contains({content, structure})) {
            return WCoefMapper::inverse_flha_mapping.at({content, structure});
        }
        LOG_ERROR("General", "Wilson coefficient with ID", content, structure, "is not supported");
    };

    static std::vector<WCoef> get_group(WGroup group) {
        switch (group) {
            case WGroup::B:
                return B_group;
            case WGroup::BPrime:
                return B_prime_group;
            case WGroup::BScalar:
                return B_scalar_group;
            case WGroup::Blnu:
                return B_lnu_group;
            case WGroup::BCLNU:
                return b_clnu_group;
            default:
                LOG_ERROR("Invalid WGroup", "get_group function couldn't find your group");
        }
    }

    static size_t n_wilsons() {
        return WCoefMapper::mapping.size();
    }

    static std::vector<std::string> get_str() {
        std::vector<std::string> _;
        for (auto&& elem : WCoefMapper::mapping) {
            _.push_back(elem.second);
        }
        return _;
    }

    static std::vector<WCoef> get_enum() {
        std::vector<WCoef> _;
        for (auto&& elem : WCoefMapper::mapping) {
            _.push_back(elem.first);
        }
        return _;
    }
private:
    static const std::vector<WCoef> B_group;
    static const std::vector<WCoef> B_prime_group;
    static const std::vector<WCoef> B_scalar_group;
    static const std::vector<WCoef> B_lnu_group;
    static const std::vector<WCoef> b_clnu_group;
    static const std::map<WCoef, std::string> mapping; 
    static const std::map<std::string, WCoef> inverse_mapping; 
    static const std::map<WCoef, std::pair<int, int>> flha_mapping; 
    static const std::map<std::pair<int, int>, WCoef> inverse_flha_mapping; 
};

/* !!!! Do not change the order of the first 4 entries !!!! */

enum class ParameterType {
    SM,
    SUSY,
    THDM,
    CUSTOM,
    FLAVOR,
    WILSON,
    DECAY,
    PASSTHROUGH,
    OBSERVABLE,
};

class ParameterTypeMapper {
public:
    static std::string str(ParameterType type) {
        return ParameterTypeMapper::mapping.at(type);
    };

    static ParameterType enum_elt(std::string type) {
        return ParameterTypeMapper::inverse_mapping.at(type);
    };

    static std::vector<std::string> get_str() {
        std::vector<std::string> _;
        for (auto&& elem : ParameterTypeMapper::mapping) {
            _.push_back(elem.second);
        }
        return _;
    }

    static std::vector<ParameterType> get_enum() {
        std::vector<ParameterType> _;
        for (auto&& elem : ParameterTypeMapper::mapping) {
            _.push_back(elem.first);
        }
        return _;
    }

private:
    static const std::map<ParameterType, std::string> mapping; 
    static const std::map<std::string, ParameterType> inverse_mapping; 
};

enum class Model {
    SM,
    SUSY,
    THDM,
    CUSTOM
};

class ModelMapper {
public:
    static std::string str(Model model) {
        return ModelMapper::mapping.at(model);
    };

    static Model enum_elt(std::string model) {
        return ModelMapper::inverse_mapping.at(model);
    };

    static std::vector<std::string> get_str() {
        std::vector<std::string> _;
        for (auto&& elem : ModelMapper::mapping) {
            _.push_back(elem.second);
        }
        return _;
    }

    static std::vector<Model> get_enum() {
        std::vector<Model> _;
        for (auto&& elem : ModelMapper::mapping) {
            _.push_back(elem.first);
        }
        return _;
    }

private:
    static const std::map<Model, std::string> mapping; 
    static const std::map<std::string, Model> inverse_mapping; 
};

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

    
    /**
     * @brief Allows for implicit conversion of a trivial LhaID to an integer 
     */
    operator long() const {
        if (this->parts.size() > 1) {
            LOG_WARN("Casting nontrivial LhaID to int discards information.");
        }
        return this->parts.at(0);
    };

    inline friend bool operator==(const LhaID& lhs, const LhaID& rhs) { return lhs.parts == rhs.parts; };
    inline friend bool operator!=(const LhaID& lhs, const LhaID& rhs) { return !(lhs == rhs); };
    inline friend bool operator<(const LhaID& lhs, const LhaID& rhs) { return (lhs.parts <=> rhs.parts) == std::weak_ordering::less; };

    friend std::ostream& operator<<(std::ostream&, const LhaID&);
};

struct ParamId {
    ParameterType type;
    std::string block;
    LhaID code;

    bool operator==(const ParamId other) const {
        return type == other.type && block == other.block && code == other.code;
    }

    bool operator<(const ParamId& other) const {
        if (type != other.type) return type < other.type;
        if (block != other.block) return block < other.block;
        return code < other.code;
    }
};

inline std::ostream& operator<<(std::ostream& os, ParamId& pid) {
    os << pid.block << ":" << pid.code;
    return os;
};

class DependenciesHelper {
public:
    static std::vector<ParamId> get_allowed_parameters(Observables id);
    static bool is_param_allowed(Observables id, ParamId pid);

private:
    static const std::map<Observables, std::vector<ParamId>> dep_lists;
};

class LhaParamsHelper {
public:
    static std::vector<std::vector<long>> get_minimal_content(const std::string& block_name);

private:
    static const std::map<std::string, std::vector<std::vector<long>>> minimal_blocks;
};


#endif // __GENERAL_H__