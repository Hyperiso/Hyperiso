#include <map>
#include <memory>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <cmath>

#include "Bs_mumu.h"
#include "Bd_mumu.h"
#include "Bu_taunu.h"

enum class ObservableEnum {
    BR_BS_MUMU,
    BR_BD_MUMU,
    BR_BU_TAUNU,
    BR_BS_MUMU_UNTAG
};

class ObservableInterface {
private:
    std::map<ObservableEnum, std::shared_ptr<Observable>> observable_map;

    void init_observable_map() {
        auto flavp = Parameters::GetInstance(3);
        double m_Bs = (*flavp)("FMASS", 531);
        double m_Bd = (*flavp)("FMASS", 511);
        double m_Bu = (*flavp)("FMASS", 521);
        observable_map = {
            {ObservableEnum::BR_BS_MUMU, std::make_shared<BR_Bs_mumu>(0, 2, m_Bs)},
            {ObservableEnum::BR_BD_MUMU, std::make_shared<BR_Bd_mumu>(0, 2, m_Bs)},
            {ObservableEnum::BR_BU_TAUNU, std::make_shared<BR_Bu_taunu>(0, 2, m_Bd)},
            {ObservableEnum::BR_BS_MUMU_UNTAG, std::make_shared<BR_Bs_mumu_untag>(0, 2, m_Bs)}
        };
    }

public:
    ObservableInterface() {
        init_observable_map();
    }

    double compute_observable(ObservableEnum obs_enum) const {
        auto it = observable_map.find(obs_enum);
        if (it != observable_map.end()) {
            return it->second->eval();
        } else {
            throw std::invalid_argument("Unknown ObservableEnum value.");
        }
    }

    double compute_variance(ObservableEnum obs_enum) const {
        auto it = observable_map.find(obs_enum);
        if (it != observable_map.end()) {
            return it->second->variance();
        } else {
            throw std::invalid_argument("Unknown ObservableEnum value.");
        }
    }

    double compute_chi2(const std::vector<ObservableEnum>& obs_list) const {
        double chi2 = 0.0;
        for (const auto& obs_enum : obs_list) {
            double value = compute_observable(obs_enum);
            double variance = compute_variance(obs_enum);
            if (variance <= 0) {
                throw std::runtime_error("Variance must be positive for chi2 calculation.");
            }
            chi2 += std::pow(value, 2) / variance;
        }
        return chi2;
    }

    ObservableEnum get_enum_from_string(const std::string& name) const {
        static const std::map<std::string, ObservableEnum> name_to_enum = {
            {"BR_BS_MUMU", ObservableEnum::BR_BS_MUMU},
            {"BR_BD_MUMU", ObservableEnum::BR_BD_MUMU},
            {"BR_BU_TAUNU", ObservableEnum::BR_BU_TAUNU},
            {"BR_BS_MUMU_UNTAG", ObservableEnum::BR_BS_MUMU}
        };
        auto it = name_to_enum.find(name);
        if (it != name_to_enum.end()) {
            return it->second;
        } else {
            throw std::invalid_argument("Invalid observable name: " + name);
        }
    }

    std::string get_string_from_enum(ObservableEnum obs_enum) const {
        static const std::map<ObservableEnum, std::string> enum_to_name = {
            {ObservableEnum::BR_BS_MUMU, "BR_BS_MUMU"},
            {ObservableEnum::BR_BD_MUMU, "BR_BD_MUMU"},
            {ObservableEnum::BR_BU_TAUNU, "BR_BU_TAUNU"},
            {ObservableEnum::BR_BS_MUMU_UNTAG, "BR_BU_TAUNU_UNTAG"}
        };
        auto it = enum_to_name.find(obs_enum);
        if (it != enum_to_name.end()) {
            return it->second;
        } else {
            throw std::invalid_argument("Unknown ObservableEnum value.");
        }
    }
};
