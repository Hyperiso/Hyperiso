#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include "Wilsonv2.h"

class CoefficientManager {
private:
    static std::map<std::string, std::unique_ptr<CoefficientManager>> instances;

    std::map<std::string, std::unique_ptr<CoefficientGroup>> coefficientGroups;

    CoefficientManager() = default;

public:
    /**
     * @brief Singleton accessor. Returns the instance of CoefficientManager for a specific model.
     * 
     * @param modelName The name of the model (e.g., "StandardModel", "THDM", "SUSY").
     * @return Pointer to the CoefficientManager instance.
     */
    static CoefficientManager* GetInstance(const std::string& modelName) {
        auto it = instances.find(modelName);
        if (it == instances.end()) {
            instances[modelName] = std::unique_ptr<CoefficientManager>(new CoefficientManager());
            return instances[modelName].get();
        }
        return it->second.get();
    }

    /**
     * @brief Register a new coefficient group under this manager.
     * 
     * @param groupName The name of the group (e.g., "BCoefficientGroup", "ScalarCoefficientGroup").
     * @param group The unique_ptr to the group to be registered.
     */
    void registerCoefficientGroup(const std::string& groupName, std::unique_ptr<CoefficientGroup> group) {
        coefficientGroups[groupName] = std::move(group);
    }

    /**
     * @brief Get a coefficient group by name.
     * 
     * @param groupName The name of the group to retrieve.
     * @return Pointer to the CoefficientGroup instance.
     */
    CoefficientGroup* getCoefficientGroup(const std::string& groupName) const {
        auto it = coefficientGroups.find(groupName);
        if (it != coefficientGroups.end()) {
            return it->second.get();
        }
        throw std::invalid_argument("CoefficientGroup not found.");
    }

    /**
     * @brief Clean up all instances.
     */
    static void Cleanup() {
        instances.clear();
    }

    /**
     * @brief Get the Wilson coefficient from a specific group.
     * 
     * @param groupName Name of the coefficient group.
     * @param coeffName Name of the coefficient.
     * @param order Order of the coefficient ("LO", "NLO", "NNLO").
     * @return The Wilson coefficient value.
     */
    std::complex<double> getCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order) const {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        auto it = group->find(coeffName);
        if (it != group->end()) {
            return it->second->get_CoefficientRunValue(order);
        }
        throw std::invalid_argument("Coefficient not found.");
    }

    /**
     * @brief Set the scale for a specific coefficient group.
     * 
     * @param groupName Name of the coefficient group.
     * @param Q The new scale to set.
     */
    void setGroupScale(const std::string& groupName, double Q) {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        group->set_Q_run(Q);
    }

    /**
     * @brief Print all coefficients of a group.
     */
    void printGroupCoefficients(const std::string& groupName) const {
        CoefficientGroup* group = getCoefficientGroup(groupName);
        std::cout << *dynamic_cast<BCoefficientGroup*>(group);
    }
};

// Initialization of the static map
std::map<std::string, std::unique_ptr<CoefficientManager>> CoefficientManager::instances;
