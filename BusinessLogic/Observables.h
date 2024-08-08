#if !defined(HYPERISO_OBSERVABLES_H)
#define HYPERISO_OBSERVABLES_H
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <iostream>

enum class Observables {
    BR_BS_MUMU,
    BR_BS_MUMU_UNTAG,
    BR_BD_MUMU,
    // BR_BU_TAUNU,
    // BR_BU_TAUNU_NP_ONLY,
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
    LAST
};

class ObservableMapper {
public:
    static ObservableMapper& getInstance() {
        static ObservableMapper instance;
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

private:

    std::unordered_map<Observables, std::string> observableToStringMap;
    ObservableMapper() {
        // Ajoutez ici les mappings initiaux
        observableToStringMap[Observables::BR_BS_MUMU] = "BR_BS_MUMU";
        observableToStringMap[Observables::BR_BS_MUMU_UNTAG] = "BR_BS_MUMU_UNTAG";
        observableToStringMap[Observables::BR_BD_MUMU] = "BR_BD_MUMU";
        observableToStringMap[Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA] = "ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA";
    }

    // Suppression des constructeurs de copie et d'assignation
    ObservableMapper(const ObservableMapper&) = delete;
    ObservableMapper& operator=(const ObservableMapper&) = delete;

    
};

#endif // HYPERISO_OBSERVABLES_H
