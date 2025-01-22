#if !defined(HYPERISO_WILSONMANAGER_H)
#define HYPERISO_WILSONMANAGER_H
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include "Wilsonv2.h"
#include "WilsonGroup.h"
#include "MemoryManager.h"
#include "QCDHelper.h"
#include <set>
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
    std::complex<double> getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};

class QSetState : public State {
public:
    QSetState(std::string order) : State(order) {this->state = StateName::QSetState;}
    void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
    void setParams(const std::string& block, int pdgCode, double value) {
        Parameters::GetInstance()->setBlockValue(block, pdgCode, value, true);
    }
    std::complex<double> getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
    std::complex<double> getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
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
    std::string model{};
    CoefficientManager() = default;
    CoefficientManager(std::string model) {this->model = model;}

public:
    // Singleton accessor
    static std::shared_ptr<CoefficientManager> GetInstance(const std::string& modelName) {
        auto it = instances.find(modelName);
        if (it == instances.end()) {
            instances[modelName] = std::shared_ptr<CoefficientManager>(new CoefficientManager(modelName));
            return instances[modelName];
        }
        return it->second;
    }

    static void initialize(const std::string& lhaFile, Model model = Model::SM, bool use_marty = false, bool is_spectrum = false, bool has_wilsons = false, bool has_obs = false) {
        MemoryManager* mm = MemoryManager::GetInstance();
        mm->init(lhaFile, model, use_marty, is_spectrum, has_wilsons, has_obs);
    }

    std::string getModel() {
        return this->model;
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
        switch(state_name) {
            case StateName::InitialState:
                groupStates[groupName] = std::make_shared<InitialState>();
                break;
            case StateName::QMatchSetState:
                groupStates[groupName] = std::make_shared<QMatchSetState>(groupStates[groupName]->getCurrentOrder());
                break;
            case StateName::MatchinSetState:
                groupStates[groupName] = std::make_shared<MatchingSetState>(groupStates[groupName]->getCurrentOrder());
                break;
            case StateName::QSetState:
                groupStates[groupName] = std::make_shared<QSetState>(groupStates[groupName]->getCurrentOrder());
                break;
            case StateName::RunSetState:
                groupStates[groupName] = std::make_shared<RunSetState>(groupStates[groupName]->getCurrentOrder());
                break;
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

    std::complex<double> getMatchingCoefficient(WilsonGroups groupName, BWilsonCoefficients coeffName, QCDOrder order) {
        return getMatchingCoefficient(GroupMapper::str(groupName), WCoefMapper::str(coeffName), OrderMapper::str(order));
    }

    std::complex<double> getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getFullMatchingCoefficient(this, groupName, coeffName, order);
    }

    std::complex<double> getFullMatchingCoefficient(WilsonGroups groupName, BWilsonCoefficients coeffName, QCDOrder order) {
        return getFullMatchingCoefficient(GroupMapper::str(groupName), WCoefMapper::str(coeffName), OrderMapper::str(order));
    }

    std::array<std::complex<double>, 3> getMatchingCoefficientSepOrders(WilsonGroups groupName, BWilsonCoefficients coeffName) {
        auto C_0 = getMatchingCoefficient(groupName, coeffName, QCDOrder::LO);
        auto C_1 = getMatchingCoefficient(groupName, coeffName, QCDOrder::NLO);
        auto C_2 = getMatchingCoefficient(groupName, coeffName, QCDOrder::NNLO);
        return {C_0, C_1, C_2};
    }

    std::vector<std::complex<double>> getAllFullMatchingCoefficients(WilsonGroups groupName, QCDOrder order) {
        std::vector<complex_t> C;
        for (auto c_id : WCoefMapper::get_group(groupName)) {
            C.emplace_back(getFullMatchingCoefficient(groupName, c_id, order));
        }
        return C;
    }

    std::vector<std::complex<double>> getAllMatchingCoefficients(WilsonGroups groupName, QCDOrder order) {
        std::vector<complex_t> C;
        for (auto c_id : WCoefMapper::get_group(groupName)) {
            C.emplace_back(getMatchingCoefficient(groupName, c_id, order));
        }
        return C;
    }

    std::complex<double> getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getRunCoefficient(this, groupName, coeffName, order);
    }

    std::complex<double> getRunCoefficient(WilsonGroups groupName, BWilsonCoefficients coeffName, QCDOrder order) {
        return getRunCoefficient(GroupMapper::str(groupName), WCoefMapper::str(coeffName), OrderMapper::str(order));
    }

    std::complex<double> getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        return ensureGroupState(groupName)->getFullRunCoefficient(this, groupName, coeffName, order);
    }

    std::complex<double> getFullRunCoefficient(WilsonGroups groupName, BWilsonCoefficients coeffName, QCDOrder order) {
        return getFullRunCoefficient(GroupMapper::str(groupName), WCoefMapper::str(coeffName), OrderMapper::str(order));
    }

    std::array<std::complex<double>, 3> getRunCoefficientSepOrders(WilsonGroups groupName, BWilsonCoefficients coeffName) {
        auto C_0 = getRunCoefficient(groupName, coeffName, QCDOrder::LO);
        auto C_1 = getRunCoefficient(groupName, coeffName, QCDOrder::NLO);
        auto C_2 = getRunCoefficient(groupName, coeffName, QCDOrder::NNLO);
        return {C_0, C_1, C_2};
    }

    std::vector<std::complex<double>> getAllFullRunCoefficients(WilsonGroups groupName, QCDOrder order) {
        std::vector<complex_t> C;
        for (auto c_id : WCoefMapper::get_group(groupName)) {
            C.emplace_back(getFullRunCoefficient(groupName, c_id, order));
        }
        return C;
    }

    std::vector<std::complex<double>> getAllRunCoefficients(WilsonGroups groupName, QCDOrder order) {
        std::vector<complex_t> C;
        for (auto c_id : WCoefMapper::get_group(groupName)) {
            C.emplace_back(getRunCoefficient(groupName, c_id, order));
        }
        return C;
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
        instances.clear();
    }

    void printGroupCoefficients(const std::string& groupName) const {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        std::cout << *dynamic_cast<BCoefficientGroup*>(group);
    }

    void update(std::string group, double Q_match=0, double Q=0) {
        Q_match = fpeq(Q_match, 0.) ? coefficientGroups[group]->get_Q_match() : Q_match; 
        Q = fpeq(Q, 0.) ? coefficientGroups[group]->get_Q_run() : Q; 
        auto order = ensureGroupState(group)->getCurrentOrder();

        this->setQMatch(group, Q_match);
        this->setMatchingCoefficient(group, order);
        this->setGroupScale(group, Q);
        this->setRunCoefficient(group, order);
    }

    static bool compatible_concat(std::shared_ptr<CoefficientManager> manager1, std::shared_ptr<CoefficientManager> manager2) {
        haveSameKeys(manager1->getGroups(), manager2->getGroups());
        for (auto& group : manager1->getGroups()) {
            if (manager1->get_state(group.first) != manager2->get_state(group.first)) {
                return false;
            }
            if (std::abs(group.second->get_Q_match() -group.second->get_Q_match()) > PRECISION) {
                return false;
            }
            if (std::abs(group.second->get_Q_run() - group.second->get_Q_run()) > PRECISION) {
                return false;
            }
        }
        return true;
    }

    static std::shared_ptr<CoefficientManager> Concat(std::shared_ptr<CoefficientManager> manager1, std::shared_ptr<CoefficientManager> manager2) {
        std::shared_ptr<CoefficientManager> newmanager = CoefficientManager::GetInstance(manager1->getModel() + "_" + manager2->getModel());
        if (!compatible_concat(manager1, manager2)) {
            LOG_ERROR("INVALID OPERATION", "impossible to concatenate", manager1->getModel(), "and", manager2->getModel());
        }
        
        for (auto& group : manager1->getGroups()) {
            std::shared_ptr<CoefficientGroup> newgroup = group.second->clone();
            for (auto& coeff : *group.second) {
                *(*newgroup)[coeff.first] += *(*manager2->getCoefficientGroup(group.first))[coeff.first];
            }
            newmanager->registerCoefficientGroup(group.first, newgroup);
            newmanager->setState(group.first, StateName::RunSetState);
        }

        return newmanager;
    }
    
    static std::shared_ptr<CoefficientManager> Builder(std::string instance, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double Q_match, double Q, std::string order) {
        std::shared_ptr<CoefficientManager> manager = CoefficientManager::GetInstance(instance);
        for (auto& group : groups) {
            manager->registerCoefficientGroup(group.first, group.second);
            manager->setQMatch(group.first, Q_match);
            manager->setMatchingCoefficient(group.first, order);
            manager->setGroupScale(group.first, Q);
            manager->setRunCoefficient(group.first, order);
        }
        return manager;
    }

    static std::shared_ptr<CoefficientManager> Builder(std::string instance, std::vector<WilsonGroups> groups, double Q_match, double Q, std::string order) {
        std::shared_ptr<CoefficientManager> manager = CoefficientManager::GetInstance(instance);
        auto group_ptrs = BuildGroupPtrs(groups);
        for (auto& group : group_ptrs) {
            manager->registerCoefficientGroup(group.first, group.second);
            manager->setQMatch(group.first, Q_match);
            manager->setMatchingCoefficient(group.first, order);
            manager->setGroupScale(group.first, Q);
            manager->setRunCoefficient(group.first, order);
        }
        return manager;
    }

private:
    static std::map<std::string, std::shared_ptr<CoefficientGroup>> BuildGroupPtrs(std::vector<WilsonGroups> group_ids) {
        std::map<std::string,std::shared_ptr<CoefficientGroup>> groups;
        for (auto g : group_ids) {
            switch (g) {
            case WilsonGroups::BCoefficients:
                groups.emplace(GroupMapper::str(g), std::make_shared<BCoefficientGroup>());
                break;
            case WilsonGroups::BPrimeCoefficients:
                groups.emplace(GroupMapper::str(g), std::make_shared<BPrimeCoefficientGroup>());
                break;
            case WilsonGroups::BScalarCoefficients:
                groups.emplace(GroupMapper::str(g), std::make_shared<BScalarCoefficientGroup>());
                break;
            }
        }
        return groups;
    }

    State* ensureGroupState(const std::string& groupName) {
        auto it = groupStates.find(groupName);
        if (it != groupStates.end()) {
            return it->second.get();
        }
        throw std::invalid_argument("State for the group not found.");
    }
};


#endif