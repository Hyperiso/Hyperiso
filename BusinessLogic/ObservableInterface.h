#include <map>
#include <memory>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <cmath>

#include "General.h"
#include "ModelEvaluator.h"

class ObservableInterface {
private:
    std::map<Observables, std::shared_ptr<Observable>> observable_map;

    void init_full_observable_map(Model model) {
        auto flavp = Parameters::GetInstance(ParameterType::FLAVOR);
        double m_Bs = (*flavp)("FMASS", 531);
        double m_Bd = (*flavp)("FMASS", 511);
        double m_Bu = (*flavp)("FMASS", 521);

        QCDOrder order = model == Model::CUSTOM ? QCDOrder::LO : QCDOrder::NNLO;
        auto B_Ks = std::make_shared<BKstarDecay>(order, Parameters::Get({ParameterType::SM, "MASS", 24}), m_Bu);

        observable_map = {
            // {Observables::BR_BS_MUMU, std::make_shared<BR_Bs_mumu>(model, order, m_Bs)},
            // {Observables::BR_BD_MUMU, std::make_shared<BR_Bd_mumu>(model, order, m_Bs)},
            // {Observables::BR_BU_TAUNU, std::make_shared<BR_Bu_taunu>(model, order, m_Bd)},
            // {Observables::BR_BS_MUMU_UNTAG, std::make_shared<BR_Bs_mumu_untag>(model, order, m_Bs)}
            {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, std::make_shared<Observable>(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, B_Ks)}
        };
    }



public:
    ObservableInterface(Model model = Model::SM) {
        init_full_observable_map(model);
    }

    ObservableInterface(const std::vector<Observables>& obs_list, Model model, QCDOrder max_order) {
        add_observables(std::move(obs_list), model, max_order);
    }

    double compute_observable(Observables obs) const {
        if (observable_map.contains(obs)) {
            return observable_map.at(obs)->eval();
        } else {
            throw std::invalid_argument("Unknown Observables value.");
        }
    }

    double compute_variance(Observables obs_enum) const {
        auto it = observable_map.find(obs_enum);
        if (it != observable_map.end()) {
            return it->second->variance();
        } else {
            throw std::invalid_argument("Unknown Observables value.");
        }
    }

    void add_observables(const std::vector<Observables>& obs_list, Model model, QCDOrder max_order) {
        if (model == Model::CUSTOM && max_order > QCDOrder::LO) {
            LOG_WARN("Only LO calculations are available on custom models, defaulting to LO");
            max_order = QCDOrder::LO;
        }

        auto B_Ks = std::make_shared<BKstarDecay>(max_order, Parameters::Get({ParameterType::SM, "MASS", 24}), Parameters::Get({ParameterType::FLAVOR, "FMASS", 521}));

        for (auto &o : obs_list) {
            switch (o) {
                // case Observables::BR_BS_MUMU:
                //     observable_map.emplace(std::make_pair(o, std::make_shared<BR_Bs_mumu>(model, max_order, p.second)));
                //     break;
                // case Observables::BR_BS_MUMU_UNTAG:
                //     observable_map.emplace(std::make_pair(o, std::make_shared<BR_Bs_mumu_untag>(model, max_order, p.second)));
                //     break;
                // case Observables::BR_BD_MUMU:
                //     observable_map.emplace(std::make_pair(o, std::make_shared<BR_Bd_mumu>(model, max_order, p.second)));
                //     break;
                // case Observables::BR_BU_TAUNU:
                //     observable_map.emplace(std::make_pair(o, std::make_shared<BR_Bu_taunu>(model, max_order, p.second)));
                //     break;
                case Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA:
                    observable_map.emplace(std::make_pair(o, std::make_shared<Observable>(o, B_Ks)));
                    break;
            }
        }
    }

    double compute_chi2() const {
        std::vector<std::shared_ptr<Observable>> obss;
        for (const auto& p : observable_map) {
            obss.emplace_back(p.second);
        }
        ModelEvaluator me (obss);
        return me.chi2();
    }
};
