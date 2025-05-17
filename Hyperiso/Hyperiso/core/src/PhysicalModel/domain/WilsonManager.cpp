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
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    this->coefficientGroups.at(groupName)->init(OrderMapper::enum_elt(order));
    if (has_bsm) {
        std::string bsm_group = groupName + bsm_suffix;
        if (!this->coefficientGroups.contains(bsm_group)) {
            throw_no_group_error(bsm_group);
        }
        this->coefficientGroups.at(bsm_group)->init(OrderMapper::enum_elt(order));
    }
}

void CoefficientManager::init_group_hadronic(const std::string& groupName, const std::string& order) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }
    this->coefficientGroups.at(groupName)->init_running_blocks(OrderMapper::enum_elt(order));
    // if (has_bsm) { //TODO : WIP with BSM
    //     std::string bsm_group = groupName + bsm_suffix;
    //     if (!this->coefficientGroups.contains(bsm_group)) {
    //         throw_no_group_error(bsm_group);
    //     }
    //     this->coefficientGroups.at(bsm_group)->init_running_blocks(OrderMapper::enum_elt(order));
    // }
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

complex_t CoefficientManager::getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, bool sm_only) {
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    complex_t c = this->coefficientGroups.at(groupName)->get_matching_coefficient(coeffName, order);

    if (has_bsm && !sm_only) {
        std::string bsm_group = groupName + bsm_suffix;
        if (!this->coefficientGroups.contains(bsm_group)) {
            throw_no_group_error(bsm_group);
        }
        c += this->coefficientGroups.at(bsm_group)->get_matching_coefficient(coeffName, order);
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
    if (!this->coefficientGroups.contains(groupName)) {
        throw_no_group_error(groupName);
    }

    complex_t c = this->coefficientGroups.at(groupName)->get_running_coefficient(coeffName, order);

    if (has_bsm && !sm_only) {
        std::string bsm_group = groupName + bsm_suffix;
        if (!this->coefficientGroups.contains(bsm_group)) {
            throw_no_group_error(bsm_group);
        }
        c += this->coefficientGroups.at(bsm_group)->get_running_coefficient(coeffName, order);
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

void CoefficientManager::post_init() {
    // TODO : Fill Wilson blocks with separate SM, BSM and SM + BSM values.
    bool marty = HAS_WILSON_API().get();
    bool SM = ModelAPI().get() == Model::SM;
    std::cout << "START INIT .................................................." << std::endl;
    if (SM) {
        if (marty) {
            for (std::string group_name : get_keys(this->coefficientGroups)) {
                WGroup group_id = GroupMapper::enum_elt(group_name);
                complete_wilson_block_from_copy(group_id, ContributionType::TOTAL, ContributionType::SM, BWilsonBasis::STANDARD);

                if (group_id == WGroup::B) {   
                    complete_wilson_block_from_copy(group_id, ContributionType::TOTAL, ContributionType::SM, BWilsonBasis::TRADITIONAL);
                }
            }
        } else {
            for (std::string group_name : get_keys(this->coefficientGroups)) {
                WGroup group_id = GroupMapper::enum_elt(group_name);
                complete_wilson_block_from_copy(group_id, ContributionType::SM, ContributionType::TOTAL, BWilsonBasis::STANDARD);

                if (group_id == WGroup::B) {   
                    complete_wilson_block_from_copy(group_id, ContributionType::SM, ContributionType::TOTAL, BWilsonBasis::TRADITIONAL);
                }
            }
        }
    } else {
        if (marty) {
            for (std::string group_name : get_keys(this->coefficientGroups)) {
                WGroup group_id = GroupMapper::enum_elt(group_name);
                complete_wilson_block_with_op(group_id, ContributionType::TOTAL, ContributionType::BSM, BWilsonBasis::STANDARD);

                if (group_id == WGroup::B) {   
                    complete_wilson_block_with_op(group_id, ContributionType::TOTAL, ContributionType::BSM, BWilsonBasis::TRADITIONAL);
                }
            }
        } else {
            for (std::string group_name : get_keys(this->coefficientGroups)) {
                // std::shared_ptr<CoefficientGroup> coefs_sm = std::make_shared<BCoefficientGroup>();
                std::string sm_name = group_name + std::string("_SM_") +  "INTER";
                std::shared_ptr<CoefficientGroup> sm_group = this->coefficientGroups[group_name]->get_sm_group();
                if (!this->coefficientGroups[group_name]) std::cerr << "null ptr\n";
                
                this->registerCoefficientGroup(sm_name, sm_group);
                this->init_group_matching(sm_name, OrderMapper::str(this->coefficientGroups[group_name]->get_order()));
                this->init_group_hadronic(sm_name, OrderMapper::str(this->coefficientGroups[group_name]->get_order()));
                WGroup group_id = GroupMapper::enum_elt(group_name);
                complete_wilson_block_with_op(group_id, ContributionType::BSM, ContributionType::TOTAL, BWilsonBasis::STANDARD);
                if (group_id == WGroup::B) {   
                    complete_wilson_block_with_op(group_id, ContributionType::BSM, ContributionType::TOTAL, BWilsonBasis::TRADITIONAL);
                }
            }
        }
    }

    // if(!SM && UseMarty) perform MARTY SM calculation to separate 0 and 1 from 2.

    // if(SM && UseMarty) copy values from 2 to 0

    // if(SM && !UseMarty) copy values from 0 to 2

    // if(!SM && !UseMarty) add 0 and 1 to build 2
}

void CoefficientManager::complete_wilson_block_with_op(WGroup group_id, ContributionType src_id, ContributionType dest_id, BWilsonBasis basis) {
    // TODO : Check whether GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis) already exists

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(group_id, ScaleType::HADRONIC) + "INTER"}},
        {ParameterType::WILSON, {GroupMapper::str(group_id, ScaleType::HADRONIC) + std::string("_SM_") + "INTER"}}
    };
    //TODO or not TODO : do not use reference if variable gonna be destroyed.
    auto func = [group_id, src_id, dest_id, basis] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        for (WCoef coef_id : WCoefMapper::get_group(group_id)) {
            for (int order=0; order<2; order++) {
                ParameterProxy pp(ParameterType::WILSON);
                complex_t coef = pp(GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis) + "INTER", WCoefMapper::flha_full(coef_id, (QCDOrder)(order + 1), src_id));
                complex_t coef_sm {};

                if (src_id != ContributionType::SM ){
                    coef_sm = pp(GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis) + std::string("SM_") + "INTER", WCoefMapper::flha_full(coef_id, (QCDOrder)(order + 1), ContributionType::SM));
                }
                // if (fpeq(coef.real(), 0.) && fpeq(coef.imag(), 0.)) continue;

                auto coef_lha_base = WCoefMapper::flha_base(coef_id);

                ParamId pid_src {
                    ParameterType::WILSON, 
                    GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), 
                    LhaID(coef_lha_base.first, coef_lha_base.second, order, (int)src_id)
                };

                ParamId pid_dest {
                    ParameterType::WILSON, 
                    GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), 
                    LhaID(coef_lha_base.first, coef_lha_base.second, order, (int)dest_id)
                };

                ParamId pid_sm {
                    ParameterType::WILSON, 
                    GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), 
                    LhaID(coef_lha_base.first, coef_lha_base.second, order, 0)
                };

                if (src_id == ContributionType::TOTAL) {
                    dep_block->store_or_assign(pid_src.code, std::make_shared<Parameter>(pid_src, coef, 0., 0.));
                    dep_block->store_or_assign(pid_dest.code, std::make_shared<Parameter>(pid_dest, coef-coef_sm, 0., 0.));
                    dep_block->store_or_assign(pid_sm.code, std::make_shared<Parameter>(pid_sm, coef_sm, 0., 0.));
                } else if (src_id == ContributionType::BSM) {
                    dep_block->store_or_assign(pid_src.code, std::make_shared<Parameter>(pid_src, coef, 0., 0.));
                    dep_block->store_or_assign(pid_dest.code, std::make_shared<Parameter>(pid_dest, coef+coef_sm, 0., 0.));
                    dep_block->store_or_assign(pid_sm.code, std::make_shared<Parameter>(pid_sm, coef_sm, 0., 0.));
                } else {
                    dep_block->store_or_assign(pid_src.code, std::make_shared<Parameter>(pid_src, coef, 0., 0.));
                    dep_block->store_or_assign(pid_dest.code, std::make_shared<Parameter>(pid_dest, coef, 0., 0.));
                }
            }
        }
    };
    WilsonParamComposer().compose_block(GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), src, func);

    std::unordered_map<ParameterType, std::vector<std::string>> src_full = {
        {ParameterType::WILSON, {GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), "WPARAM_RUN_SM"}}
    };
    // std::cout << this->coefficientGroups.at(GroupMapper::str(group_id)) << std::endl;
    this->coefficientGroups.at(GroupMapper::str(group_id))->init_full_running_block(src_full, basis, false, {src_id, dest_id}); //TODO : We need to add full again ?
    // this->coefficientGroups.at(GroupMapper::str(group_id))->init_full_running_block(src_full, basis, false, dest_id); //TODO : We need to add full again ?
}

void CoefficientManager::complete_wilson_block_from_copy(WGroup group_id, ContributionType src_id, ContributionType dest_id, BWilsonBasis basis) {
    // TODO : Check whether GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis) already exists

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(group_id, ScaleType::HADRONIC) + "INTER"}}
    };
    //TODO or not TODO : do not use reference if variable gonna be destroyed.
    auto func = [group_id, src_id, dest_id, basis] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        for (WCoef coef_id : WCoefMapper::get_group(group_id)) {
            for (int order=0; order<2; order++) {
                ParameterProxy pp(ParameterType::WILSON);
                complex_t coef = pp(GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis) + "INTER", WCoefMapper::flha_full(coef_id, (QCDOrder)(order + 1), src_id));
                
                if (fpeq(coef.real(), 0.) && fpeq(coef.imag(), 0.)) continue;

                auto coef_lha_base = WCoefMapper::flha_base(coef_id);
                ParamId pid_src {
                    ParameterType::WILSON, 
                    GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), 
                    LhaID(coef_lha_base.first, coef_lha_base.second, order, (int)src_id)
                };

                ParamId pid_dest {
                    ParameterType::WILSON, 
                    GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), 
                    LhaID(coef_lha_base.first, coef_lha_base.second, order, (int)dest_id)
                };

                dep_block->store_or_assign(pid_src.code, std::make_shared<Parameter>(pid_src, coef, 0., 0.));
                dep_block->store_or_assign(pid_dest.code, std::make_shared<Parameter>(pid_dest, coef, 0., 0.));
            }
        }
    };
    WilsonParamComposer().compose_block(GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), src, func);

    std::unordered_map<ParameterType, std::vector<std::string>> src_full = {
        {ParameterType::WILSON, {GroupMapper::str(group_id, ScaleType::HADRONIC, false, basis), "WPARAM_RUN_SM"}}
    };
    // std::cout << this->coefficientGroups.at(GroupMapper::str(group_id)) << std::endl;
    this->coefficientGroups.at(GroupMapper::str(group_id))->init_full_running_block(src_full, basis, false, {src_id, dest_id}); //TODO : We need to add full again ?
    // this->coefficientGroups.at(GroupMapper::str(group_id))->init_full_running_block(src_full, basis, false, dest_id); //TODO : We need to add full again ?
}

std::shared_ptr<CoefficientManager> CoefficientManager::Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order) {
    if (groups.empty()) {
        LOG_WARN("(CoefficientManager) No coefficient groups provided.");
        return std::make_shared<CoefficientManager>();
    }

    WilsonParameterHelper().init(2);
    auto manager = std::make_shared<CoefficientManager>();
    manager->has_bsm = model == ModelMapper::str(Model::THDM) || model == ModelMapper::str(Model::SUSY);
    manager->bsm_suffix = manager->has_bsm ? "_" + model : "";
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
    }
    manager->post_init();
    LOG_INFO("(CoefficientManager) Manager successfully initialized");
    return manager;
}

void CoefficientManager::set_hadronic_scale(double mu_h) {
    ScaleSetter(ScaleType::HADRONIC).set(mu_h);
}

void CoefficientManager::set_matching_scale(double mu_W) {
    ScaleSetter(ScaleType::MATCHING).set(mu_W);
}
