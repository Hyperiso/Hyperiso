#include "WilsonManager.h"

void CoefficientManager::initialize(const std::string& lhaFile, Model model, bool use_marty, bool is_spectrum, bool has_wilsons, bool has_obs) {
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

std::string CoefficientManager::getModel() {
    return ModelMapper::str(ModelAPI().get());
}

void CoefficientManager::init_group_matching(const std::string& groupName, const std::string& order) {
    std::cout << "trying to init" << std::endl;
    this->coefficientGroups.at(groupName)->init(OrderMapper::enum_elt(order));
    std::cout << "trying to init " << groupName<< std::endl;
    if (has_bsm) {
        this->coefficientGroups.at(groupName + bsm_suffix)->init(OrderMapper::enum_elt(order));
    }
}

void CoefficientManager::init_group_hadronic(const std::string& groupName, const std::string& order) {
    this->coefficientGroups.at(groupName)->init_running_block(OrderMapper::enum_elt(order));
    if (has_bsm) {
        this->coefficientGroups.at(groupName + bsm_suffix)->init_running_block(OrderMapper::enum_elt(order));
    }
}

void CoefficientManager::switchbasis(const std::string& groupName) {
    this->coefficientGroups.at(groupName)->switch_basis();
    if (has_bsm) {
        this->coefficientGroups.at(groupName + bsm_suffix)->switch_basis();
    }
}

complex_t CoefficientManager::getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only) {
    complex_t c = this->coefficientGroups.at(groupName)->get_matching_coefficient(coeffName, order);
    if (has_bsm && !sm_only) {
        c += this->coefficientGroups.at(groupName + bsm_suffix)->get_matching_coefficient(coeffName, order);
    }
    return c;
}

complex_t CoefficientManager::getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only) {
    double fact = wilson_p("WPARAM_MATCH_SM", 1) / (4 * PI);
    int max_order = static_cast<int>(OrderMapper::enum_elt(order));
    complex_t c {0};
    for (size_t o = 1; o <= max_order; o++) {
        c += this->getMatchingCoefficient(groupName, coeffName, OrderMapper::str(static_cast<QCDOrder>(o)), sm_only) * std::pow(fact, o - 1);
    }
    return c;
}

complex_t CoefficientManager::getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only) {
    complex_t c = this->coefficientGroups.at(groupName)->get_running_coefficient(coeffName, order);
    if (has_bsm && !sm_only) {
        c += this->coefficientGroups.at(groupName + bsm_suffix)->get_running_coefficient(coeffName, order);
    }
    return c;
}

complex_t CoefficientManager::getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only) {
    double fact = wilson_p("WPARAM_RUN_SM", 1) / (4 * PI);
    int max_order = static_cast<int>(OrderMapper::enum_elt(order));
    complex_t c {0};
    for (size_t o = 1; o <= max_order; o++) {
        c += this->getRunCoefficient(groupName, coeffName, OrderMapper::str(static_cast<QCDOrder>(o)), sm_only) * std::pow(fact, o - 1);
    }
    return c;
}

void CoefficientManager::registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group) {
    coefficientGroups[groupName] = group;
}

std::shared_ptr<CoefficientGroup> CoefficientManager::getCoefficientGroup(const std::string& groupName) const {
    auto it = coefficientGroups.find(groupName);
    if (it != coefficientGroups.end()) {
        return it->second;
    }
    throw std::invalid_argument("CoefficientGroup not found.");
}

std::map<std::string, std::shared_ptr<CoefficientGroup>> CoefficientManager::getGroups() {
    return this->coefficientGroups;
}

void CoefficientManager::printGroupCoefficients(const std::string& groupName) const {
    std::shared_ptr<CoefficientGroup> group = getCoefficientGroup(groupName);
    std::cout << group;
}

void CoefficientManager::update(std::string group, double mu_W, double mu_h) {
    this->set_matching_scale(mu_W);
    this->set_hadronic_scale(mu_h);
}

CoefficientManager CoefficientManager::Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order) {
    if (groups.empty()) {
        LOG_WARN("(CoefficientManager) No coefficient groups provided.");
        return CoefficientManager();
    }

    WilsonParameterHelper().init(2);
    CoefficientManager manager;
    manager.has_bsm = model == ModelMapper::str(Model::THDM) || model == ModelMapper::str(Model::SUSY);
    manager.bsm_suffix = manager.has_bsm ? "_" + model : "";
    for (auto& group : groups) {
        LOG_INFO("(CoefficientManager) Registering coefficient group", group.first);
        manager.registerCoefficientGroup(group.first, group.second);
    }
    LOG_INFO("(CoefficientManager) Setting matching scale");
    manager.set_matching_scale(mu_W);
    LOG_INFO("(CoefficientManager) Setting hadronic scale");
    manager.set_hadronic_scale(mu_h);
    for (auto& group: groups) {
        if (manager.has_bsm && group.first.ends_with(manager.bsm_suffix)) continue;
        LOG_INFO("(CoefficientManager) Initializing group matching", group.first);
        manager.init_group_matching(group.first, order);
        LOG_INFO("(CoefficientManager) Initializing group hadronic", group.first);
        manager.init_group_hadronic(group.first, order);
    }
    return manager;
}

void CoefficientManager::set_hadronic_scale(double mu_h) {
    this->hss.set(mu_h);
}

void CoefficientManager::set_matching_scale(double mu_W) {
    this->mss.set(mu_W);
}
