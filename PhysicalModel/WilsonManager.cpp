#include "WilsonManager.h"

InitialState::InitialState() {
    this->state = StateName::InitialState;
}

InitialState::~InitialState() {

}

void InitialState::setQMatch(CoefficientManager* manager, const std::string& groupName, double Q_match) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_match(Q_match);
        manager->setState(groupName, std::make_shared<QMatchSetState>(OrderMapper::str(this->currentOrder)));
}

void MatchingSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        group->set_Q_run(Q);
        manager->setState(groupName, std::make_shared<QSetState>(OrderMapper::str(this->currentOrder)));
}


void QMatchSetState::setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        QCDOrder newOrder = order == "LO" ? QCDOrder::LO :
                                    order == "NLO" ? QCDOrder::NLO :
                                    QCDOrder::NNLO;

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
        manager->setState(groupName, std::make_shared<MatchingSetState>(order));
}

void MatchingSetState::setMatchingCoefficient(CoefficientManager* manager, const std::string& groupName, const std::string& order) {
        CoefficientGroup* group = manager->getCoefficientGroup(groupName);
        QCDOrder newOrder = order == "LO" ? QCDOrder::LO :
                                    order == "NLO" ? QCDOrder::NLO :
                                    QCDOrder::NNLO;
        // Wilson_parameters::GetInstance()->SetMuW(manager->getCoefficientGroup(groupName)->get_Q_match());
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

double MatchingSetState::getAlphaS(CoefficientManager* manager, const std::string& groupName) {
    return QCDHelper::alpha_s(manager->getCoefficientGroup(groupName)->get_Q_match());
}

double QMatchSetState::getAlphaS(CoefficientManager* manager, const std::string& groupName) {
    return QCDHelper::alpha_s(manager->getCoefficientGroup(groupName)->get_Q_match());
}

void RunSetState::setGroupScale(CoefficientManager* manager, const std::string& groupName, double Q) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    group->set_Q_run(Q);
    if (OrderMapper::str(currentOrder) == "LO") {
        group->set_base_1_LO();
    } else if (OrderMapper::str(currentOrder) == "NLO") {
        group->set_base_1_NLO();
    } else if (OrderMapper::str(currentOrder) == "NNLO") {
        group->set_base_1_NNLO();
    }
}

void RunSetState::setQMatch(CoefficientManager *manager,
                            const std::string &groupName,
                            double Q_match) {
    CoefficientGroup* group = manager->getCoefficientGroup(groupName);
    group->set_Q_match(Q_match);
    manager->setState(groupName, std::make_shared<QMatchSetState>(OrderMapper::str(this->currentOrder)));
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
    group->switch_base();
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

std::map<std::string, std::shared_ptr<CoefficientManager>> CoefficientManager::instances;