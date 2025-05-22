#include "WilsonManager.h"

void CoefficientManager::throw_no_group_error(const std::string &groupName) const {
    std::stringstream ss;
    ss << "Coefficient group " << groupName << " not found in manager. Existing groups are:\n";
    for (const auto& group : this->coefficientGroups) {
        ss << "\t- " << group.first << "\n";
    }
    LOG_ERROR("KeyError", ss.str());
}

void CoefficientManager::initialize(const std::string &lhaFile,
                                    Model model,
                                    bool use_marty,
                                    bool is_spectrum,
                                    bool has_wilsons,
                                    bool has_obs)
{
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
    int only_total = 0;
    switch (OrderMapper::enum_elt(order))
    {
    case QCDOrder::NNLO:
        init_specific_order_group_matching(groupName, OrderMapper::str(QCDOrder::NNLO), only_total++);
        [[fallthrough]];
        
    case QCDOrder::NLO:
        init_specific_order_group_matching(groupName, OrderMapper::str(QCDOrder::NLO), only_total++);
        [[fallthrough]];

    case QCDOrder::LO:
        init_specific_order_group_matching(groupName, OrderMapper::str(QCDOrder::LO), only_total++);
        break;

    default:
        break;
    }
}

void CoefficientManager::init_specific_order_group_matching(const std::string& groupName, const std::string& order, bool only_total) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    bool marty = HAS_WILSON_API().get();
    bool SM = ModelAPI().get() == Model::SM;
    
    if (!only_total) {
        LOG_INFO("Initializing group", groupName, "at", order);
        this->coefficientGroups.at(groupName)->init(OrderMapper::enum_elt(order));
        if (!SM) {
            LOG_INFO("Computing SM contribution");
            std::string sm_group = groupName + "_SM";
            std::shared_ptr<CoefficientGroup> sm_group_ptr = this->coefficientGroups[groupName]->get_sm_group();
            if (!sm_group_ptr)
                LOG_ERROR("LogicError", "No SM group found for " + groupName);
            
            this->registerCoefficientGroup(sm_group, sm_group_ptr);
            sm_group_ptr->init(OrderMapper::enum_elt(order));
        }
    }


    QCDOrder enum_order = OrderMapper::enum_elt(order);
    std::string storage_block = this->coefficientGroups.at(groupName)->get_matching_storage_block();
    WilsonParamComposer composer;

    for (auto& coeff : *this->coefficientGroups.at(groupName)) {
        if (SM) {
            ContributionType c_type = marty ? ContributionType::SM : ContributionType::TOTAL;
            ParamId pid {
                storage_block,
                WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff.second->get_base_name()), enum_order, c_type)
            };

            composer.compose_parameter(
                pid, 
                std::unordered_set<ParamId> {{storage_block, coeff.second->id(enum_order, coeff.second->get_type())}}, 
                [&] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
                    dep_param->set_expected(src.at({ParameterType::WILSON, storage_block, coeff.second->id(enum_order, coeff.second->get_type())})->get_val());
                }
            );
        } else {

            ContributionType c_type = marty ? ContributionType::BSM : ContributionType::TOTAL;
            ParamId pid_dest {
                this->coefficientGroups.at(groupName)->get_matching_storage_block(),
                WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff.second->get_base_name()), enum_order, c_type)
            };

            ParamId pid_src {
                this->coefficientGroups.at(groupName)->get_matching_storage_block(),
                WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff.second->get_base_name()), enum_order, ContributionType::SM)
            };

            std::unordered_set<ParamId> sources {pid_src, {storage_block, coeff.second->id(enum_order, coeff.second->get_type())}};

            composer.compose_parameter(
                pid_dest, 
                sources, 
                [&] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
                    double src_val = ParameterProxy(ParameterType::WILSON)(storage_block, coeff.second->id(enum_order, coeff.second->get_type()));
                    double sm_val = ParameterProxy(ParameterType::WILSON)(pid_src.block, pid_src.code);
                    dep_param->set_expected(marty ? src_val - sm_val : src_val + sm_val);
                }
            );
        }
    }

    LOG_INFO("CoefficientManager", "Initialized group", groupName, "at", order);
}

void CoefficientManager::fill_sources_for_group(const std::string & groupName, const std::string& order, std::unordered_map<ParameterType, std::vector<std::string>>& src, int id) {
    switch (OrderMapper::enum_elt(order))
    {
    case QCDOrder::NNLO:
        for (const auto& [key, vec] : this->coefficientGroups[groupName]->get_sources(QCDOrder::NNLO, id)) {
            auto& targetVec = src[key];
            targetVec.insert(targetVec.end(), vec.begin(), vec.end());
        }
        [[fallthrough]];
        
    case QCDOrder::NLO:
        for (const auto& [key, vec] : this->coefficientGroups[groupName]->get_sources(QCDOrder::NLO, id)) {
            auto& targetVec = src[key];
            targetVec.insert(targetVec.end(), vec.begin(), vec.end());
        }
        [[fallthrough]];

    case QCDOrder::LO:
        for (const auto& [key, vec] : this->coefficientGroups[groupName]->get_sources(QCDOrder::LO, id)) {
            auto& targetVec = src[key];
            targetVec.insert(targetVec.end(), vec.begin(), vec.end());
        }
        break;

    default:
        break;
    }
}

std::pair<WCoef, std::pair<QCDOrder, ContributionType>> lha_wilson_deserialize(LhaID id) {
    auto parts = id.get_parts();
    auto w_id = std::make_pair<int, int>(parts[0], parts[1]);

    if (!WCoefMapper::inverse_flha_mapping().contains(w_id)) {
        LOG_ERROR("ValueError", "bad lha id for wilson conversion");
    }
    WCoef coef = WCoefMapper::inverse_flha_mapping().at(w_id);
    QCDOrder order = parts[2] ? ((parts[2] -1) ? QCDOrder::NNLO : QCDOrder::NLO) : QCDOrder::LO;
    ContributionType part = parts[3] ? parts[3] -1 ? ContributionType::TOTAL : ContributionType::BSM : ContributionType::SM;

    std::pair<WCoef, std::pair<QCDOrder, ContributionType>> ret;
    ret = {coef, {order, part}};

    return ret;
}

//TODO : manage id correctly, use base ENUM.
void CoefficientManager::init_group_hadronic(const std::string& groupName, const std::string& order, int id) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    std::unordered_map<ParameterType, std::vector<std::string>> src = {};
    fill_sources_for_group(groupName, order, src, id);

    QCDOrder ord = OrderMapper::enum_elt(order);

    std::map<QCDOrder, std::function<std::unordered_map<WCoef, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&, const std::unordered_map<std::string, std::shared_ptr<Block>>&)>> funcs = {
        {QCDOrder::LO, this->coefficientGroups[groupName]->get_func(QCDOrder::LO, id)},
        {QCDOrder::NLO, this->coefficientGroups[groupName]->get_func(QCDOrder::NLO, id)},
        {QCDOrder::NNLO, this->coefficientGroups[groupName]->get_func(QCDOrder::NNLO, id)}
    };

    std::string matching_block_name = this->coefficientGroups[groupName]->get_matching_storage_block();

    auto func = [matching_block_name, ord, funcs, groupName] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        std::map<LhaID, std::shared_ptr<Parameter>> matching_coeff = src.at(matching_block_name)->getItems();
        std::cout << src.at(matching_block_name) << std::endl;
        std::unordered_map<ContributionType, std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>> matching_map;
        for (auto& coef : matching_coeff) {
            std::pair<WCoef, std::pair<QCDOrder, ContributionType>> c = lha_wilson_deserialize(coef.first);

            const WCoef& wcoef = c.first;
            const QCDOrder& order = c.second.first;
            const ContributionType& contrib = c.second.second;
            matching_map[contrib][order][wcoef] = coef.second->get_val();

        }
        std::unordered_map<ContributionType, std::unordered_map<QCDOrder,std::unordered_map<WCoef, scalar_t>>> res;
        for (auto contri : {ContributionType::SM, ContributionType::BSM, ContributionType::TOTAL}) {
            switch (ord)
                {
                case QCDOrder::NNLO:
                    res[contri][QCDOrder::NNLO] = funcs.at(QCDOrder::NNLO)(matching_map[contri], src);
                    [[fallthrough]];
                    
                case QCDOrder::NLO:
                    res[contri][QCDOrder::NLO] = funcs.at(QCDOrder::NLO)(matching_map[contri], src);
                    [[fallthrough]];

                case QCDOrder::LO:
                    res[contri][QCDOrder::LO] = funcs.at(QCDOrder::LO)(matching_map[contri], src);
                    break;

                default:
                    break;
                }
        }
        for (auto& [c_type, order_map] : res) { // Iterate over the contributions
            for (auto& [order, coef_map] : order_map) { // Iterate over the orders
                for (auto& [coef_id, coef_val] : coef_map) { // Iterate over the coefficients
                    LhaID coef_lha = WCoefMapper::flha_full(coef_id, order, c_type);
                    ParamId pid {
                        ParameterType::WILSON, 
                        GroupMapper::str(GroupMapper::enum_elt(groupName), ScaleType::HADRONIC, false, BWilsonBasis::STANDARD), 
                        coef_lha
                    };
                    dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, coef_val, 0., (int)c_type));
                }
            }
        }
    };

    WilsonParamComposer().compose_block(GroupMapper::str(GroupMapper::enum_elt(groupName), ScaleType::HADRONIC, false, id ? BWilsonBasis::TRADITIONAL : BWilsonBasis::STANDARD), src, func);
}

void CoefficientManager::switchbasis(const std::string& groupName) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    this->coefficientGroups.at(groupName)->switch_basis();
    if (has_bsm) {
        std::string bsm_group = groupName + bsm_suffix;
        if (!this->coefficientGroups.contains(bsm_group)) {
            throw_no_group_error(bsm_group);
        }
        this->coefficientGroups.at(bsm_group)->switch_basis();
    }
}

complex_t CoefficientManager::getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    complex_t c = this->coefficientGroups.at(groupName)->get_matching_coefficient(coeffName, order, cont_type);

    return c;
}

complex_t CoefficientManager::getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type) {
    double fact = wilson_p("WPARAM_MATCH_SM", 1) / (4 * PI);
    int max_order = static_cast<int>(OrderMapper::enum_elt(order));
    complex_t c {0};
    for (size_t o = 1; o <= max_order; o++) {
        c += this->getMatchingCoefficient(groupName, coeffName, OrderMapper::str(static_cast<QCDOrder>(o)), cont_type) * std::pow(fact, o - 1);
    }
    return c;
}

complex_t CoefficientManager::getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    complex_t c = this->coefficientGroups.at(groupName)->get_running_coefficient(coeffName, order, cont_type);
    return c;
}

complex_t CoefficientManager::getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type) {
    double fact = wilson_p("WPARAM_RUN_SM", 1) / (4 * PI);
    int max_order = static_cast<int>(OrderMapper::enum_elt(order));
    complex_t c {0};
    for (size_t o = 1; o <= max_order; o++) {
        c += this->getRunCoefficient(groupName, coeffName, OrderMapper::str(static_cast<QCDOrder>(o)), cont_type) * std::pow(fact, o - 1);
    }
    return c;
}

void CoefficientManager::registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group) {
    coefficientGroups[groupName] = group;
}

std::shared_ptr<CoefficientGroup> CoefficientManager::getCoefficientGroup(const std::string& groupName) const {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    return this->coefficientGroups.at(groupName);
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

std::shared_ptr<CoefficientManager> CoefficientManager::Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order) {
    if (groups.empty()) {
        LOG_WARN("(CoefficientManager) No coefficient groups provided.");
        return std::make_shared<CoefficientManager>();
    }

    WilsonParameterHelper().init(2);
    auto manager = std::make_shared<CoefficientManager>();
    bool is_bsm = model != ModelMapper::str(Model::SM);
    for (auto& group : groups) {
        LOG_INFO("(CoefficientManager) Registering coefficient group", group.first);
        manager->registerCoefficientGroup(group.first, group.second);
    }
    LOG_INFO("(CoefficientManager) Setting matching scale");
    manager->set_matching_scale(mu_W);
    LOG_INFO("(CoefficientManager) Setting hadronic scale");
    manager->set_hadronic_scale(mu_h);
    for (auto& group: groups) {
        if (manager->has_bsm && group.first.ends_with(manager->bsm_suffix)) continue;
        LOG_INFO("(CoefficientManager) Initializing group matching", group.first);
        manager->init_group_matching(group.first, order);
        LOG_INFO("(CoefficientManager) Initializing group hadronic", group.first);
        manager->init_group_hadronic(group.first, order);
        if (group.second->is_double_basis()){
            manager->init_group_hadronic(group.first, order, 1); //TODO : base 2}
        }
    }
    LOG_INFO("(CoefficientManager) Manager successfully initialized");
    return manager;
}

void CoefficientManager::set_hadronic_scale(double mu_h) {
    ScaleSetter(ScaleType::HADRONIC).set(mu_h);
}

void CoefficientManager::set_matching_scale(double mu_W) {
    ScaleSetter(ScaleType::MATCHING).set(mu_W);
}
