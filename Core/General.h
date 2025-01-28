#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <string>
#include <stdexcept>
#include <unordered_map>
#include <iostream>
#include <map>
#include <vector>
#include "Logger.h"

enum class Observables {
    BR_BS_MUMU,
    BR_BS_MUMU_UNTAG,
    BR_BD_MUMU,
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
    BR_B_XS_GAMMA
};

class ObservableMapper {
public:
    static std::string str(Observables obs) {
        return ObservableMapper::mapping.at(obs);
    };

    static Observables enum_elt(std::string name) {
        return ObservableMapper::inverse_mapping.at(name);
    };

private:
    static const std::map<Observables, std::string> mapping; 
    static const std::map<std::string, Observables> inverse_mapping; 
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

private:
    static const std::map<QCDOrder, std::string> mapping; 
    static const std::map<std::string, QCDOrder> inverse_mapping; 
};

enum class WCoef {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2
};

enum class WGroup {
    B, 
    BPrime, 
    BScalar
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

    static std::string flha(WCoef coef) {
        return WCoefMapper::flha_mapping.at(coef);
    };

    static WCoef from_flha(std::string flha_id) {
        if (WCoefMapper::inverse_flha_mapping.contains(flha_id)) {
            return WCoefMapper::inverse_flha_mapping.at(flha_id);
        }
        LOG_ERROR("General", "Wilson coefficient with ID", flha_id, "is not supported");
    };

    static std::vector<WCoef> get_group(WGroup group) {
        switch (group) {
            case WGroup::B:
                return B_group;
            case WGroup::BPrime:
                return B_prime_group;
            case WGroup::BScalar:
                return B_scalar_group;
        }
    }

    static size_t n_wilsons() {
        return WCoefMapper::mapping.size();
    }

private:
    static const std::vector<WCoef> B_group;
    static const std::vector<WCoef> B_prime_group;
    static const std::vector<WCoef> B_scalar_group;
    static const std::map<WCoef, std::string> mapping; 
    static const std::map<std::string, WCoef> inverse_mapping; 
    static const std::map<WCoef, std::string> flha_mapping; 
    static const std::map<std::string, WCoef> inverse_flha_mapping; 
};

/* !!!! Do not change the order of the first 4 entries !!!! */

enum class ParameterType {
    SM,
    SUSY,
    THDM,
    CUSTOM,
    FLAVOR,
    WILSON,
    FF
};

class ParameterTypeMapper {
public:
    static std::string str(ParameterType type) {
        return ParameterTypeMapper::mapping.at(type);
    };

    static ParameterType enum_elt(std::string type) {
        return ParameterTypeMapper::inverse_mapping.at(type);
    };

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

private:
    static const std::map<Model, std::string> mapping; 
    static const std::map<std::string, Model> inverse_mapping; 
};

struct ParamId {
    ParameterType type;
    std::string block;
    int code;

    bool operator==(const ParamId other) const {
        return type == other.type && block == other.block && code == other.code;
    }

    bool operator<(const ParamId& other) const {
        if (type != other.type) return type < other.type;
        if (block != other.block) return block < other.block;
        return code < other.code;
    }
};

#endif // __GENERAL_H__