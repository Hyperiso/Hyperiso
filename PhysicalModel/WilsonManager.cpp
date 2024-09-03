#include "WilsonManager.h"

InitialState::InitialState() {
    this->state = "InitialState";
}

InitialState::~InitialState() {

}

void InitialState::setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_match(Q_match);
        std::cout << "Q match set." << std::endl;
        manager->setState(groupName, std::make_unique<QMatchSetState>(this->EnumToString(this->currentOrder)));
}

void MatchingSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_run(Q);
        std::cout << "Scale set to " << Q << " for group: " << groupName << std::endl;
        manager->setState(groupName, std::make_unique<QSetState>(this->EnumToString(this->currentOrder)));
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
                it.second->LO_calculation();
                it.second->NLO_calculation();
            } else if (order == "NNLO") {
                it.second->LO_calculation();
                it.second->NLO_calculation();
                it.second->NNLO_calculation();
            }
        }

        currentOrder = newOrder;
        std::cout << "Matching coefficients set for group: " << groupName << ", order: " << order << std::endl;
        manager->setState(groupName, std::make_unique<MatchingSetState>(order));
}

void QSetState::setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
    std::cout << "order is LO ? " << (currentOrder == CoefficientOrder::LO) << std::endl;
    std::cout << "order is NLO ? " << (currentOrder == CoefficientOrder::NLO) << std::endl;
    std::cout << "order is nNLO ? " << (currentOrder == CoefficientOrder::NNLO) << std::endl;
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
    manager->setState(groupName, std::make_unique<RunSetState>(this->EnumToString(this->currentOrder)));
}

std::complex<double> MatchingSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getMatching(coeffName, order);
}

std::complex<double> MatchingSetState::getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullMatching(coeffName, order);
}

void RunSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    group->set_Q_run(Q);
    std::cout << "Scale set to " << Q << " for group: " << groupName << std::endl;
    if (this->EnumToString(currentOrder) == "LO") {
        group->set_base_1_LO();
    } else if (this->EnumToString(currentOrder) == "NLO") {
        group->set_base_1_NLO();
    } else if (this->EnumToString(currentOrder) == "NNLO") {
        group->set_base_1_NNLO();
    }
}

std::complex<double> RunSetState::getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullMatching(coeffName, order);
}

std::complex<double> QSetState::getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullMatching(coeffName, order);
}

void RunSetState::switchbasis(CoefficientManager* manager, const std::string& groupName) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->switch_base();
}
std::complex<double> RunSetState::getFullRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullRun(coeffName, order);
}

std::complex<double> RunSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getMatching(coeffName, order);

}

std::complex<double> QSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getMatching(coeffName, order);

}

std::complex<double> RunSetState::getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Run coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getRun(coeffName, order);
}

// Initialization of the static map
std::map<std::string, std::unique_ptr<CoefficientManager>> CoefficientManager::instances;