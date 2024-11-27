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
    BR_BU_TAUNU,
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA
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
    LO,
    NLO,
    NNLO,
    NONE,
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

enum class BWilsonCoefficients {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2
};

enum class WilsonGroups {
    BCoefficients, 
    BPrimeCoefficients, 
    BScalarCoefficients
};

enum class BWilsonBasis {
    STANDARD, 
    TRADITIONAL
};

class GroupMapper {
public:
    static std::string str(WilsonGroups group) {
        return GroupMapper::mapping.at(group);
    };

    static WilsonGroups enum_elt(std::string group) {
        return GroupMapper::inverse_mapping.at(group);
    };

private:
    static const std::map<WilsonGroups, std::string> mapping; 
    static const std::map<std::string, WilsonGroups> inverse_mapping; 
};

class WCoefMapper {

public:
    static std::string str(BWilsonCoefficients coef) {
        return WCoefMapper::mapping.at(coef);
    };

    static BWilsonCoefficients enum_elt(std::string coef) {
        return WCoefMapper::inverse_mapping.at(coef);
    };

    static std::string flha(BWilsonCoefficients coef) {
        return WCoefMapper::flha_mapping.at(coef);
    };

    static BWilsonCoefficients from_flha(std::string flha_id) {
        if (WCoefMapper::inverse_flha_mapping.contains(flha_id)) {
            return WCoefMapper::inverse_flha_mapping.at(flha_id);
        }
        LOG_ERROR("General", "Wilson coefficient with ID", flha_id, "is not supported");
    };

    static std::vector<BWilsonCoefficients> get_group(WilsonGroups group) {
        switch (group) {
            case WilsonGroups::BCoefficients:
                return B_group;
            case WilsonGroups::BPrimeCoefficients:
                return B_prime_group;
            case WilsonGroups::BScalarCoefficients:
                return B_scalar_group;
        }
    }

    static size_t n_wilsons() {
        return WCoefMapper::mapping.size();
    }

private:
    static const std::vector<BWilsonCoefficients> B_group;
    static const std::vector<BWilsonCoefficients> B_prime_group;
    static const std::vector<BWilsonCoefficients> B_scalar_group;
    static const std::map<BWilsonCoefficients, std::string> mapping; 
    static const std::map<std::string, BWilsonCoefficients> inverse_mapping; 
    static const std::map<BWilsonCoefficients, std::string> flha_mapping; 
    static const std::map<std::string, BWilsonCoefficients> inverse_flha_mapping; 
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

#endif // __GENERAL_H__