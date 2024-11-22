#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <string>
#include <stdexcept>
#include <unordered_map>
#include <iostream>
#include <map>
#include <vector>

enum class Observables {
    FIRST,
    BR_BS_MUMU,
    BR_BS_MUMU_UNTAG,
    BR_BD_MUMU,
    BR_BU_TAUNU,
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
    LAST
};

class ObservableMapper {
    static ObservableMapper* instance;
public:

    static ObservableMapper* GetInstance() {
        if (!instance) {
            instance = new ObservableMapper();
        }
        return instance;
    }

    void addMapping(Observables observable, const std::string& name) {
        observableToStringMap[observable] = name;
    }

    std::string getString(Observables observable) const {
        auto it = observableToStringMap.find(observable);
        if (it != observableToStringMap.end()) {
            return it->second;
        } else {
            throw std::runtime_error("Observable not found in map");
        }
    }

    Observables getObservable(const std::string& name) const {
        for (const auto& pair : observableToStringMap) {
            if (pair.second == name) {
                return pair.first;
            }
        }
        throw std::runtime_error("String not found in map");
    }

    std::unordered_map<Observables, std::string> get_map() {
        return this->observableToStringMap;
    }

private:

    std::unordered_map<Observables, std::string> observableToStringMap;
    ObservableMapper() {
        observableToStringMap[Observables::BR_BS_MUMU] = "BR_Bsmumu";
        observableToStringMap[Observables::BR_BS_MUMU_UNTAG] = "BRuntag_Bsmumu";
        observableToStringMap[Observables::BR_BD_MUMU] = "BR_Bdmumu";
        observableToStringMap[Observables::BR_BU_TAUNU] = "BR_Bu_Taunu";
        observableToStringMap[Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA] = "ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA";
    }

    ObservableMapper(const ObservableMapper&) = delete;
    ObservableMapper& operator=(const ObservableMapper&) = delete;
};

enum class CoefficientOrder {
    NONE,
    LO,
    NLO,
    NNLO
};

class OrderMapper {

public:
    static std::string str(CoefficientOrder order) {
        return OrderMapper::mapping.at(order);
    };

    static CoefficientOrder enum_elt(std::string order) {
        return OrderMapper::inverse_mapping.at(order);
    };

private:
    static const std::map<CoefficientOrder, std::string> mapping; 
    static const std::map<std::string, CoefficientOrder> inverse_mapping; 
};

enum class WilsonCoefficientList {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2
};

enum class WilsonGroups {
    BCoefficients, 
    BPrimeCoefficients, 
    BScalarCoefficients
};

class GroupMapper {
public:
    static std::string str(WilsonGroups order) {
        return GroupMapper::mapping.at(order);
    };

    static WilsonGroups enum_elt(std::string order) {
        return GroupMapper::inverse_mapping.at(order);
    };

private:
    static const std::map<WilsonGroups, std::string> mapping; 
    static const std::map<std::string, WilsonGroups> inverse_mapping; 
};

class WCoefMapper {

public:
    static std::string str(WilsonCoefficientList order) {
        return WCoefMapper::mapping.at(order);
    };

    static WilsonCoefficientList enum_elt(std::string order) {
        return WCoefMapper::inverse_mapping.at(order);
    };

    static std::string flha(WilsonCoefficientList order) {
        return WCoefMapper::flha_mapping.at(order);
    };

    static std::vector<WilsonCoefficientList> get_group(WilsonGroups group) {
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
    static const std::vector<WilsonCoefficientList> B_group;
    static const std::vector<WilsonCoefficientList> B_prime_group;
    static const std::vector<WilsonCoefficientList> B_scalar_group;
    static const std::map<WilsonCoefficientList, std::string> mapping; 
    static const std::map<std::string, WilsonCoefficientList> inverse_mapping; 
    static const std::map<WilsonCoefficientList, std::string> flha_mapping; 
};


#endif // __GENERAL_H__