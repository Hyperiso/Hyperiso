#if !defined(HYPERISO_WILSONMANAGER_H)
#define HYPERISO_WILSONMANAGER_H
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include "Wilsonv2.h"

class CoefficientManager; // Forward declaration

enum class CoefficientOrder {
    NONE,
    LO,
    NLO,
    NNLO
};

class State {
protected:
    CoefficientOrder currentOrder = CoefficientOrder::NONE;

public:
    virtual ~State() = default;

    virtual void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
        throw std::runtime_error("Invalid state: Cannot set group scale in current state.");
    }

    virtual void setQMatch(CoefficientManager* manager) {
        throw std::runtime_error("Invalid state: Cannot set Q match in current state.");
    }

    virtual void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot set matching coefficients in current state.");
    }

    virtual void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot set run coefficients in current state.");
    }

    virtual void getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get matching coefficients in current state.");
    }

    virtual void getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
        throw std::runtime_error("Invalid state: Cannot get run coefficients in current state.");
    }

    bool isOrderCalculated(const std::string& order) {
        if (order == "LO" && currentOrder >= CoefficientOrder::LO) {
            return true;
        }
        if (order == "NLO" && currentOrder >= CoefficientOrder::NLO) {
            return true;
        }
        if (order == "NNLO" && currentOrder >= CoefficientOrder::NNLO) {
            return true;
        }
        return false;
    }
};

// class QStetState;
// class QMatchSetState;
// class MatchingSetState;
// class RunSetState;

class InitialState : public State {
public:
    void setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) override;

    InitialState();
    ~InitialState();
};


class QSetState : public State {
public:
    void setQMatch(CoefficientManager* manager) override;
};


class QMatchSetState : public State {
public:
    void setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;
};


class MatchingSetState : public State {
public:
    void setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) override;

    void getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};

class RunSetState : public State {
public:
    void getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;

    void getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) override;
};



class CoefficientManager {
private:
    static std::map<std::string, std::unique_ptr<CoefficientManager>> instances;
    std::map<std::string, std::unique_ptr<CoefficientGroup>> coefficientGroups;
    std::unique_ptr<State> state;

    CoefficientManager() : state(std::make_unique<InitialState>()) {}

public:
    // Singleton accessor
    static CoefficientManager* GetInstance(const std::string& modelName) {
        auto it = instances.find(modelName);
        if (it == instances.end()) {
            instances[modelName] = std::unique_ptr<CoefficientManager>(new CoefficientManager());
            return instances[modelName].get();
        }
        return it->second.get();
    }

    inline void setState(std::unique_ptr<State> newState) {
        state = std::move(newState);
    }

    inline void setGroupScale(const std::string& groupName, double Q) {
        state->setGroupScale(this, groupName, Q);
    }

    void setQMatch() {
        state->setQMatch(this);
    }

    inline void setMatchingCoefficient(const std::string& groupName, const std::string& order) {
        state->setMatchingCoefficient(this, groupName, order);
    }

    inline void setRunCoefficient(const std::string& groupName, const std::string& order) {
        state->setRunCoefficient(this, groupName, order);
    }

    inline void getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        state->getMatchingCoefficient(this, groupName, coeffName, order);
    }

    inline void getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) {
        state->getRunCoefficient(this, groupName, coeffName, order);
    }

    // Register a new coefficient group under this manager
    inline void registerCoefficientGroup(const std::string& groupName, std::unique_ptr<CoefficientGroup> group) {
        coefficientGroups[groupName] = std::move(group);
    }

    // Get a coefficient group by name
    CoefficientGroup* getCoefficientGroup(const std::string& groupName) const {
        auto it = coefficientGroups.find(groupName);
        if (it != coefficientGroups.end()) {
            return it->second.get();
        }
        throw std::invalid_argument("CoefficientGroup not found.");
    }

    // Clean up all instances
    static void Cleanup() {
        instances.clear();
    }

    // Print all coefficients of a group
    inline void printGroupCoefficients(const std::string& groupName) const {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        std::cout << *dynamic_cast<BCoefficientGroup*>(group);
    }
};




// class CoefficientManager {
// private:
//     static std::map<std::string, std::unique_ptr<CoefficientManager>> instances;

//     std::map<std::string, std::unique_ptr<CoefficientGroup>> coefficientGroups;


//     CoefficientManager() = default;

// public:
//     /**
//      * @brief Singleton accessor. Returns the instance of CoefficientManager for a specific model.
//      * 
//      * @param modelName The name of the model (e.g., "StandardModel", "THDM", "SUSY").
//      * @return Pointer to the CoefficientManager instance.
//      */
//     static CoefficientManager* GetInstance(const std::string& modelName) {
//         auto it = instances.find(modelName);
//         if (it == instances.end()) {
//             instances[modelName] = std::unique_ptr<CoefficientManager>(new CoefficientManager());
//             return instances[modelName].get();
//         }
//         return it->second.get();
//     }

//     /**
//      * @brief Register a new coefficient group under this manager.
//      * 
//      * @param groupName The name of the group (e.g., "BCoefficientGroup", "ScalarCoefficientGroup").
//      * @param group The unique_ptr to the group to be registered.
//      */
//     void registerCoefficientGroup(const std::string& groupName, std::unique_ptr<CoefficientGroup> group) {
//         coefficientGroups[groupName] = std::move(group);
//     }

//     /**
//      * @brief Get a coefficient group by name.
//      * 
//      * @param groupName The name of the group to retrieve.
//      * @return Pointer to the CoefficientGroup instance.
//      */
//     CoefficientGroup* getCoefficientGroup(const std::string& groupName) const {
//         auto it = coefficientGroups.find(groupName);
//         if (it != coefficientGroups.end()) {
//             return it->second.get();
//         }
//         throw std::invalid_argument("CoefficientGroup not found.");
//     }

//     /**
//      * @brief Clean up all instances.
//      */
//     static void Cleanup() {
//         instances.clear();
//     }

//     /**
//      * @brief Get the Wilson coefficient from a specific group.
//      * 
//      * @param groupName Name of the coefficient group.
//      * @param coeffName Name of the coefficient.
//      * @param order Order of the coefficient ("LO", "NLO", "NNLO").
//      * @return The Wilson coefficient value.
//      */
//     std::complex<double> getCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) const {
//         CoefficientGroup* group = getCoefficientGroup(groupName);
//         auto it = group->find(coeffName);
//         if (it != group->end()) {
//             return it->second->get_CoefficientRunValue(order);
//         }
//         throw std::invalid_argument("Coefficient not found.");
//     }

//     /**
//      * @brief Set the matching Wilson coefficients from a specific group.
//      * 
//      * @param groupName Name of the coefficient group.
//      * @param order Order of the coefficient ("LO", "NLO", "NNLO").
//      * @return None.
//      */
//     void setMatchingCoefficient(const std::string& groupName, const std::string& order) const {
//         CoefficientGroup* group = getCoefficientGroup(groupName);
//         for (auto& it : *group) {
//             if (order == "LO")
//                 it.second->LO_calculation();
//             if (order == "NLO")
//                 it.second->NLO_calculation();
//             if (order == "NNLO")
//                 it.second->NNLO_calculation();
//         }
//     }

//     /**
//      * @brief Set the Run Wilson coefficient from a specific group.
//      * 
//      * @param groupName Name of the coefficient group.
//      * @param order Order of the coefficient ("LO", "NLO", "NNLO").
//      * @return None.
//      */
//     std::complex<double> setRunCoefficient(const std::string& groupName, const std::string& order) const {
//         CoefficientGroup* group = getCoefficientGroup(groupName);
//         if (order == "LO")
//             group->set_base_1_LO();
//         if (order == "NLO")
//             group->set_base_1_NLO();
//         if (order == "NNLO")
//             group->set_base_1_NNLO();
//     }
//     /**
//      * @brief Set the scale for a specific coefficient group.
//      * 
//      * @param groupName Name of the coefficient group.
//      * @param Q The new scale to set.
//      */
//     void setGroupScale(const std::string& groupName, double Q) {
//         CoefficientGroup* group = getCoefficientGroup(groupName);
//         group->set_Q_run(Q);
//     }

//     /**
//      * @brief Print all coefficients of a group.
//      */
//     void printGroupCoefficients(const std::string& groupName) const {
//         CoefficientGroup* group = getCoefficientGroup(groupName);
//         std::cout << *dynamic_cast<BCoefficientGroup*>(group);
//     }
// };

// // Initialization of the static map
// std::map<std::string, std::unique_ptr<CoefficientManager>> CoefficientManager::instances;

#endif