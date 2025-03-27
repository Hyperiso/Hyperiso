#include "WilsonManager.h"

InitialState::InitialState() {
    this->state = StateName::InitialState;
}

InitialState::~InitialState() {

}

void InitialState::setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        std::cout << "eheh" << std::endl;
        group->set_Q_match(Q_match);
        std::cout << "eheh2" << std::endl;
        manager->setState(groupName, std::make_shared<QMatchSetState>(OrderMapper::str(this->currentOrder)));
}

void MatchingSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_run(Q);
        manager->setState(groupName, std::make_shared<QSetState>(OrderMapper::str(this->currentOrder)));
}

void MatchingSetState::setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_match(Q_match);
        manager->setState(groupName, std::make_shared<QMatchSetState>(OrderMapper::str(this->currentOrder)));
}

void QMatchSetState::setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        QCDOrder newOrder = OrderMapper::enum_elt(order);

        if (newOrder < currentOrder) {
            throw std::runtime_error("Cannot set matching coefficient: Lower or same order already calculated.");
        }

        for (auto& it : *group) {
            if (order == "LO") {
                LOG_INFO("Performing LO calculation of Wilson coefficient");
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
        manager->setState(groupName, std::make_shared<MatchingSetState>(order));
}

void MatchingSetState::setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        QCDOrder newOrder = OrderMapper::enum_elt(order);

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
        manager->setState(groupName, std::make_shared<MatchingSetState>(order));
}

void QSetState::setRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    if (order == "LO") {
        group->set_base_1_LO();
    } else if (order == "NLO") {
        group->set_base_1_LO();
        group->set_base_1_NLO();
    } else if (order == "NNLO") {
        group->set_base_1_LO();
        group->set_base_1_NLO();
        group->set_base_1_NNLO();
    }

    manager->setState(groupName, std::make_unique<RunSetState>(OrderMapper::str(this->currentOrder)));
}

complex_t MatchingSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getMatching(coeffName, order);
}

complex_t MatchingSetState::getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullMatching(coeffName, order);
}

double MatchingSetState::getAlphaS(CoefficientManager* manager, const std::string& groupName) {
    return QCDHelper::alpha_s(manager->getCoefficientGroup(groupName)->get_Q_match());
}

double QMatchSetState::getAlphaS(CoefficientManager* manager, const std::string& groupName) {
    return QCDHelper::alpha_s(manager->getCoefficientGroup(groupName)->get_Q_match());
}

void RunSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    group->set_Q_run(Q);
    
    switch(currentOrder) {
        case QCDOrder::LO:
            group->set_base_1_LO();
            break;
        case QCDOrder::NLO:
            group->set_base_1_NLO();
            break;
        case QCDOrder::NNLO:
            group->set_base_1_NNLO();
            break;
    }
}

void RunSetState::setQMatch(CoefficientManager *manager,
                            const std::string &groupName,
                            double Q_match) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    group->set_Q_match(Q_match);
    manager->setState(groupName, std::make_shared<QMatchSetState>(OrderMapper::str(this->currentOrder)));
}

complex_t RunSetState::getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullMatching(coeffName, order);
}

complex_t QSetState::getFullMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullMatching(coeffName, order);
}

void RunSetState::switchbasis(CoefficientManager* manager, const std::string& groupName) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    group->switch_base();
}
complex_t RunSetState::getFullRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getfullRun(coeffName, order);
}

complex_t RunSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getMatching(coeffName, order);

}

complex_t QSetState::getMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Matching coefficient of the requested order has not been set.");
    }
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getMatching(coeffName, order);

}

complex_t RunSetState::getRunCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& coeffName, const std::string& order) {
    
    if (!isOrderCalculated(order)) {
        throw std::runtime_error("Run coefficient of the requested order has not been set.");
    }

    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    return group->getRun(coeffName, order);
}

std::shared_ptr<CoefficientManager> CoefficientManager::instance;