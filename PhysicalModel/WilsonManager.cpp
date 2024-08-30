#include "WilsonManager.h"

InitialState::InitialState() {

}

InitialState::~InitialState() {

}

void InitialState::setQMatch(CoefficientManager* manager) {
        std::cout << "Q match set." << std::endl;
        manager->setState(std::make_unique<QMatchSetState>());
}

void MatchingSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_run(Q);
        std::cout << "Scale set to " << Q << " for group: " << groupName << std::endl;
        manager->setState(std::make_unique<QSetState>());
}



void QMatchSetState::setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        CoefficientOrder newOrder = order == "LO" ? CoefficientOrder::LO :
                                    order == "NLO" ? CoefficientOrder::NLO :
                                    CoefficientOrder::NNLO;

        if (newOrder <= currentOrder) {
            throw std::runtime_error("Cannot set matching coefficient: Lower or same order already calculated.");
        }

        for (auto& it : *group) {
            if (order == "LO") {
                it.second->LO_calculation();
            } else if (order == "NLO") {
                it.second->NLO_calculation();
            } else if (order == "NNLO") {
                it.second->NNLO_calculation();
            }
        }

        currentOrder = newOrder;
        std::cout << "Matching coefficients set for group: " << groupName << ", order: " << order << std::endl;
        manager->setState(std::make_unique<MatchingSetState>());
}

void QSetState::setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    if (order == "LO") {
        group->set_base_1_LO();
    } else if (order == "NLO") {
        group->set_base_1_NLO();
    } else if (order == "NNLO") {
        group->set_base_1_NNLO();
    }

    std::cout << "Run coefficients set for group: " << groupName << ", order: " << order << std::endl;
    manager->setState(std::make_unique<RunSetState>());
}

void MatchingSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    auto it = group->find(coeffName);
    if (it != group->end()) {
        std::cout << "Matching coefficient: " << it->second->get_CoefficientRunValue(order) << std::endl;
    } else {
        throw std::invalid_argument("Matching Coefficient not found.");
    }
}

void RunSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    auto it = group->find(coeffName);
    if (it != group->end()) {
        std::cout << "Matching coefficient: " << it->second->get_CoefficientRunValue(order) << std::endl;
    } else {
        throw std::invalid_argument("Matching Coefficient not found.");
    }
}

void RunSetState::getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Run coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    auto it = group->find(coeffName);
    if (it != group->end()) {
        std::cout << "Run coefficient: " << it->second->get_CoefficientMatchingValue(order) << std::endl;
    } else {
        throw std::invalid_argument("Run Coefficient not found.");
    }
}

// Initialization of the static map
std::map<std::string, std::unique_ptr<CoefficientManager>> CoefficientManager::instances;