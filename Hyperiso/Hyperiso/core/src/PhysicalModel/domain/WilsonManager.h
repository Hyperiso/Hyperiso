#if !defined(HYPERISO_WILSONMANAGER_H)
#define HYPERISO_WILSONMANAGER_H
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <set>
#include "Wilson.h"
#include "WilsonGroup.h"
#include "MemoryManager.h"
#include "QCDHelper.h"
#include "Utils.h"
#include "HyperisoMaster.h"
#define PRECISION 1e-5

inline bool haveSameKeys(const std::map<std::string, std::shared_ptr<CoefficientGroup>>& map1, const std::map<std::string, std::shared_ptr<CoefficientGroup>>& map2) {
    std::set<std::string> keys1;
    std::set<std::string> keys2;

    for (const auto& pair : map1) {
        keys1.insert(pair.first);
    }

    for (const auto& pair : map2) {
        keys2.insert(pair.first);
    }

    return keys1 == keys2;
}

enum class StateName {
    InitialState, QMatchSetState, MatchinSetState, QSetState, RunSetState
};

class CoefficientManager; // Forward declaration

class State {
protected:
    QCDOrder currentOrder = QCDOrder::NONE;
    StateName state{};
public:
    virtual ~State() = default;

    State() = default;
    State(std::string order) {
        currentOrder = OrderMapper::enum_elt(order);
    }

    virtual void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
        throw std::runtime_error("Invalid state: Cannot set group scale in current state.");
    }

    virtual void setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) {
        throw std::runtime_error("Invalid state: Cannot set Q match in current state.");
    }
    virtual void setParams(const std::string& block, int pdgCode, double value) {
        throw std::runtime_error("Invalid state: Cannot set params in current state.");
    }

    virtual void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot set matching coefficients in current state.");
    }

    virtual void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot set run coefficients in current state.");
    }

    virtual void switchbasis(CoefficientManager* manager, const std::string& groupName) {
        throw std::runtime_error("Invalid state: Cannot switch basis if not already calculated");
    }

    virtual complex_t getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get matching coefficients in current state.");
    }

    virtual complex_t getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get full matching coefficients in current state.");
    }

    virtual complex_t getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get run coefficients in current state.");
    }

    virtual complex_t getFullRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get full run coefficients in current state.");
    }

    virtual double getAlphaS(CoefficientManager* manager, const std::string& groupName)  {
        throw std::runtime_error("Invalid state: Cannot get alpha_s in current state.");
    }
    
    bool isOrderCalculated(const std::string& order) {
        if (order == "LO" && currentOrder >= QCDOrder::LO) {
            return true;
        }
        if (order == "NLO" && currentOrder >= QCDOrder::NLO) {
            return true;
        }
        if (order == "NNLO" && currentOrder >= QCDOrder::NNLO) {
            return true;
        }
        return false;
    }
    
    std::string getCurrentOrder() {
        return OrderMapper::str(this->currentOrder);
    }

    StateName get_state() {
        return this->state;
    }
};


class InitialState : public State {
public:
    void setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) override;

    InitialState();
    ~InitialState();

};


class QMatchSetState : public State {
public:
    QMatchSetState(std::string order) : State(order) { this->state = StateName::QMatchSetState;}
    void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    double getAlphaS(CoefficientManager* manager, const std::string& groupName);
};


class MatchingSetState : public State {
public:
    MatchingSetState(std::string order) : State(order) {this->state = StateName::MatchinSetState;}
    void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) override;
    // void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) override;
    void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    double getAlphaS(CoefficientManager* manager, const std::string& groupName);
    complex_t getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    complex_t getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};

class QSetState : public State {
public:
    QSetState(std::string order) : State(order) {this->state = StateName::QSetState;}
    void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    complex_t getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    complex_t getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};

class RunSetState : public State {
public:
    RunSetState(std::string order) : State(order) {this->state = StateName::RunSetState;}
    void switchbasis(CoefficientManager* manager, const std::string& groupName);
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) override;
    void setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) override;
    complex_t getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    complex_t getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    complex_t getFullRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    complex_t getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};



class CoefficientManager {
private:
    static std::shared_ptr<CoefficientManager> instance;
    std::map<std::string, std::shared_ptr<CoefficientGroup>> coefficientGroups;
    std::map<std::string, std::shared_ptr<State>> groupStates;

    bool has_bsm;
    std::string bsm_suffix;

    CoefficientManager() = default;

public:
    static std::shared_ptr<CoefficientManager> GetInstance() {
        if (!instance) {
            instance = std::shared_ptr<CoefficientManager>(new CoefficientManager());
        }
        return instance;
    }

    static void initialize(const std::string& lhaFile, Model model = Model::SM, bool use_marty = false, bool is_spectrum = false, bool has_wilsons = false, bool has_obs = false) {
        MemoryManager* mm = MemoryManager::GetInstance();
        HyperisoMaster hyp = HyperisoMaster(); //TODO bad coupling
        Config config;
        config.flags.at(ExternalFlag::IS_LHA_SPECTRUM) = is_spectrum;
        config.flags.at(ExternalFlag::USE_MARTY) = use_marty;
        config.flags.at(ExternalFlag::HAS_WILSON_INPUT) = has_wilsons;
        config.flags.at(ExternalFlag::HAS_TH_OBSERVABLE_INPUT) = has_obs;

        config.model = model;
    
        hyp.init(lhaFile, config);
    }

    std::string getModel() {
        return ModelMapper::str(ModelAPI().get());
    }

    double get_params(const std::string& block, int pdgCode) {
        return (*Parameters::GetInstance())(block, pdgCode);
    }

    StateName get_state(const std::string& groupName) {
        return ensureGroupState(groupName)->get_state();
    }

    void setState(const std::string& groupName, std::shared_ptr<State> newState) {
        groupStates[groupName] = std::move(newState);
    }

    void setState(const std::string& groupName, StateName state_name) {
        std::shared_ptr<State> new_state;
        [[maybe_unused]] auto order = groupStates[groupName]->getCurrentOrder();
        switch(state_name) {
            case StateName::InitialState:
                new_state = std::make_shared<InitialState>();
                break;
            case StateName::QMatchSetState:
                new_state = std::make_shared<QMatchSetState>(order);
                break;
            case StateName::MatchinSetState:
                new_state = std::make_shared<MatchingSetState>(order);
                break;
            case StateName::QSetState:
                new_state = std::make_shared<QSetState>(order);
                break;
            case StateName::RunSetState:
                new_state = std::make_shared<RunSetState>(order);
                break;
        }

        setState(groupName, new_state);
    }

    void setGroupScale(const std::string& groupName, double Q) {
        ensureGroupState(groupName)->setGroupScale(this, groupName, Q);
        if (has_bsm) {
            ensureGroupState(groupName + bsm_suffix)->setGroupScale(this, groupName + bsm_suffix, Q);
        }
    }

    void setQMatch(const std::string& groupName, double Q_match) {
        std::cout << "here 1" << std::endl;
        ensureGroupState(groupName)->setQMatch(this, groupName, Q_match);
        std::cout << "here 2" << std::endl;
        if (has_bsm) {
            std::cout << "here 3" << std::endl;
            ensureGroupState(groupName + bsm_suffix)->setQMatch(this, groupName + bsm_suffix, Q_match);
        }
    }

    void setParams(const std::string& groupName, const std::string& block, int pdgCode, double value) {
        ensureGroupState(groupName)->setParams(block, pdgCode, value);
    }

    void setMatchingCoefficient(const std::string& groupName, const std::string& order) {
        ensureGroupState(groupName)->setMatchingCoefficient(this, groupName, order);
        if (has_bsm) {
            ensureGroupState(groupName + bsm_suffix)->setMatchingCoefficient(this, groupName + bsm_suffix, order);
        }
    }

    void setRunCoefficient(const std::string& groupName, const std::string& order) {
        ensureGroupState(groupName)->setRunCoefficient(this, groupName, order);
        if (has_bsm) {
            ensureGroupState(groupName + bsm_suffix)->setRunCoefficient(this, groupName + bsm_suffix, order);
        }
    }

    void switchbasis(const std::string& groupName) {
        ensureGroupState(groupName)->switchbasis(this, groupName);
        if (has_bsm) {
            ensureGroupState(groupName + bsm_suffix)->switchbasis(this, groupName + bsm_suffix);
        }
    }
    
    complex_t getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false) {
        complex_t c = ensureGroupState(groupName)->getMatchingCoefficient(this, groupName, coeffName, order);
        if (has_bsm && !sm_only) {
            c += ensureGroupState(groupName + bsm_suffix)->getMatchingCoefficient(this, groupName + bsm_suffix, coeffName, order);
        }
        return c;
    }

    complex_t getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false) {
        complex_t c = ensureGroupState(groupName)->getFullMatchingCoefficient(this, groupName, coeffName, order);
        if (has_bsm && !sm_only) {
            c += ensureGroupState(groupName + bsm_suffix)->getFullMatchingCoefficient(this, groupName + bsm_suffix, coeffName, order);
        }
        return c;
    }

    complex_t getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false) {
        complex_t c = ensureGroupState(groupName)->getRunCoefficient(this, groupName, coeffName, order);
        if (has_bsm && !sm_only) {
            c += ensureGroupState(groupName + bsm_suffix)->getRunCoefficient(this, groupName + bsm_suffix, coeffName, order);
        }
        return c;
    }

    complex_t getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only=false) {
        complex_t c = ensureGroupState(groupName)->getFullRunCoefficient(this, groupName, coeffName, order);
        if (has_bsm && !sm_only) {
            c += ensureGroupState(groupName + bsm_suffix)->getFullRunCoefficient(this, groupName + bsm_suffix, coeffName, order);
        }
        return c;
    }

    double getAlphaS(const std::string& groupName) {
        return ensureGroupState(groupName)->getAlphaS(this, groupName);
    }
    
    void registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group) {
        coefficientGroups[groupName] = group;
        groupStates[groupName] = std::make_shared<InitialState>();
    }

    CoefficientGroup* getCoefficientGroup(const std::string& groupName) const {
        auto it = coefficientGroups.find(groupName);
        if (it != coefficientGroups.end()) {
            return it->second.get();
        }
        throw std::invalid_argument("CoefficientGroup not found.");
    }

    std::map<std::string, std::shared_ptr<CoefficientGroup>> getGroups() {
        return this->coefficientGroups;
    }

    static void Cleanup() {
        instance.reset();
    }

    void printGroupCoefficients(const std::string& groupName) const {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        std::cout << *dynamic_cast<BCoefficientGroup*>(group);
    }

    void update(std::string group, double Q_match=0, double Q=0) {
        auto apply_update = [this, Q_match, Q] (std::string wg) { 
            auto order = ensureGroupState(wg)->getCurrentOrder();
            this->setQMatch(wg, Q_match);
            this->setMatchingCoefficient(wg, order);
            this->setGroupScale(wg, Q);
            this->setRunCoefficient(wg, order);
        };

        apply_update(group);
        if (has_bsm) {
            apply_update(group + bsm_suffix);
        }
    }
    
    static std::shared_ptr<CoefficientManager> Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double Q_match, double Q, std::string order) {
        std::shared_ptr<CoefficientManager> manager = CoefficientManager::GetInstance();
        manager->has_bsm = model == ModelMapper::str(Model::THDM) || model == ModelMapper::str(Model::SUSY);
        manager->bsm_suffix = manager->has_bsm ? "_" + model : "";
        std::cout << "fuck 1" << std::endl;
        for (auto& group : groups) {
            LOG_DEBUG("(CoefficientManager) Registering coefficient group", group.first);
            manager->registerCoefficientGroup(group.first, group.second);
        }
        std::cout << "fuck 2" << std::endl;
        for (auto& group: groups) {
            if (manager->has_bsm && group.first.ends_with(manager->bsm_suffix)) continue;

            manager->setQMatch(group.first, Q_match);
            std::cout << "fuck 2.5" << std::endl;
            manager->setMatchingCoefficient(group.first, order);
            std::cout << "fuck 2.8" << std::endl;
            manager->setGroupScale(group.first, Q);
            std::cout << "fuck 2.9" << std::endl;
            manager->setRunCoefficient(group.first, order);
        }
        std::cout << "fuck 3" << std::endl;

        return manager;
    }

private:
    State* ensureGroupState(const std::string& groupName) {
        auto it = groupStates.find(groupName);
        if (it != groupStates.end()) {
            return it->second.get();
        }
        throw std::invalid_argument("State for the group not found.");
    }
};


#endif