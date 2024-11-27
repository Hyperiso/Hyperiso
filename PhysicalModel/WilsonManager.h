#if !defined(HYPERISO_WILSONMANAGER_H)
#define HYPERISO_WILSONMANAGER_H
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include "Wilsonv2.h"
#include "MemoryManager.h"

class CoefficientManager; // Forward declaration

class State {
protected:
    QCDOrder currentOrder = QCDOrder::NONE;
    std::string state{};
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

    virtual std::complex<double> getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get matching coefficients in current state.");
    }

    virtual std::complex<double> getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get full matching coefficients in current state.");
    }

    virtual std::complex<double> getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get run coefficients in current state.");
    }

    virtual std::complex<double> getFullRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
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

    std::string get_state() {
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
    QMatchSetState(std::string order) : State(order) { this->state = "QMatchSetState";}
    void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    double getAlphaS(CoefficientManager* manager, const std::string& groupName);
};


class MatchingSetState : public State {
public:
    MatchingSetState(std::string order) : State(order) {this->state = "MatchingSetState";}
    void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) override;
    // void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    double getAlphaS(CoefficientManager* manager, const std::string& groupName);
    std::complex<double> getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};

class QSetState : public State {
public:
    QSetState(std::string order) : State(order) {this->state = "QSetState";}
    void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    std::complex<double> getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};

class RunSetState : public State {
public:
    RunSetState(std::string order) : State(order) {this->state = "RunSetState";}
    void switchbasis(CoefficientManager* manager, const std::string& groupName);
    void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) override;
    std::complex<double> getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getFullRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};



class CoefficientManager {
private:
    static std::map<std::string, std::shared_ptr<CoefficientManager>> instances;
    std::map<std::string, std::shared_ptr<CoefficientGroup>> coefficientGroups;
    std::map<std::string, std::shared_ptr<State>> groupStates;

    CoefficientManager() = default;

public:
    // Singleton accessor
    static CoefficientManager* GetInstance(const std::string& modelName) {
        auto it = instances.find(modelName);
        if (it == instances.end()) {
            instances[modelName] = std::shared_ptr<CoefficientManager>(new CoefficientManager());
            return instances[modelName].get();
        }
        return it->second.get();
    }

    static void initialize(const std::string& lhaFile, Model model) {
        MemoryManager* mm = MemoryManager::GetInstance();
        mm->init(lhaFile, model);
    }

    double get_params(const std::string& block, int pdgCode) {
        return (*Parameters::GetInstance())(block, pdgCode);
    }
    std::string get_state(const std::string& groupName) {
        return ensureGroupState(groupName)->get_state();
    }

    void setState(const std::string& groupName, std::shared_ptr<State> newState) {
        groupStates[groupName] = std::move(newState);
    }

    void setState(const std::string& groupName, std::string state_name) {
        if (state_name == "matchingstate"){
            groupStates[groupName] = std::make_shared<MatchingSetState>(groupStates[groupName]->getCurrentOrder());
        }
    }
    void setGroupScale(const std::string& groupName, double Q) {
        ensureGroupState(groupName)->setGroupScale(this, groupName, Q);
    }

    void setQMatch(const std::string& groupName, double Q_match) {
        ensureGroupState(groupName)->setQMatch(this, groupName, Q_match);
    }

    void setParams(const std::string& groupName, const std::string& block, int pdgCode, double value) {
        ensureGroupState(groupName)->setParams(block, pdgCode, value);
    }
    void setMatchingCoefficient(const std::string& groupName, const std::string& order) {
        ensureGroupState(groupName)->setMatchingCoefficient(this, groupName, order);
    }

    void setRunCoefficient(const std::string& groupName, const std::string& order) {
        ensureGroupState(groupName)->setRunCoefficient(this, groupName, order);
    }

    void switchbasis(const std::string& groupName) {
        ensureGroupState(groupName)->switchbasis(this, groupName);
    }
    std::complex<double> getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getMatchingCoefficient(this, groupName, coeffName, order);
    }

    std::complex<double> getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getFullMatchingCoefficient(this, groupName, coeffName, order);
    }

    std::complex<double> getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getRunCoefficient(this, groupName, coeffName, order);
    }

    std::complex<double> getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getFullRunCoefficient(this, groupName, coeffName, order);
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

    static void Cleanup() {
        instances.clear();
    }

    void printGroupCoefficients(const std::string& groupName) const {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        std::cout << *dynamic_cast<BCoefficientGroup*>(group);
    }

    static CoefficientManager* Builder(std::string instance, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double Q_match, double Q, std::string order) {
        CoefficientManager *manager = CoefficientManager::GetInstance(instance);
        for (auto& group : groups) {
            manager->registerCoefficientGroup(group.first, group.second);
            manager->setQMatch(group.first, Q_match);
            manager->setMatchingCoefficient(group.first, order);
            manager->setGroupScale(group.first, Q);
            manager->setRunCoefficient(group.first, order);
        }
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