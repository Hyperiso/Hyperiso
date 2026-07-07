#include "WilsonBuilder.h"

#include "GroupDefinition.h"
#include "WilsonCoefficientRegistry.h"
#include "CoefficientGroupBuilder.h"

WilsonBuilder::WilsonBuilder(WilsonBuildConfig config) {
    this->build(config);
}

static std::shared_ptr<CoefficientRegistry> make_registry() {
    auto reg = std::make_shared<CoefficientRegistry>();
    register_B(*reg);
    register_BPrime(*reg);
    register_BScalar(*reg);
    register_CC_bc(*reg);
    register_CC_bu(*reg);
    register_CC_cs(*reg);
    register_CC_cd(*reg);
    register_CC_su(*reg);
    register_CC_du(*reg);
    register_K(*reg);
    register_MesonMixing(*reg);
    return reg;
}

WilsonBuilder::WilsonBuilder(std::shared_ptr<CoefficientManager> manager) : cm(manager) {}

void WilsonBuilder::build(WilsonBuildConfig config) {

    std::shared_ptr<IBlockComposer> iblock_c = std::make_shared<WilsonParamComposer>();
    
    wilson_param_helpers[Model::SM] = std::make_shared<WilsonParameterHelper>(iblock_c);

    for (const auto& elem : config.groups) {
        wilson_param_helpers[Model::SM]->init(2, elem);

    }
    
    Model model = ModelAPI().get();

    if (model == Model::THDM) {
        wilson_param_helpers[model] = std::make_shared<THDMParameterHelper>(iblock_c);
        for (const auto& elem : config.groups) {
            wilson_param_helpers[model]->init(2, elem);

        }
    } else if (model == Model::SUSY) {
        wilson_param_helpers[model] = std::make_shared<SUSYParameterHelper>(iblock_c);
        for (const auto& elem : config.groups) {
            wilson_param_helpers[model]->init(2, elem);
        }
    }

    

    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy = std::make_shared<ParameterProxy>(ParameterType::WILSON);
    std::shared_ptr<IParameterProxy<std::string, LhaID>> sm_proxy = std::make_shared<ParameterProxy>(ParameterType::SM);
    std::shared_ptr<ICoreAPI<bool>> use_marty = std::make_shared<UseMarty>();
    std::shared_ptr<ICoreAPI<bool>> has_wilson = std::make_shared<HasWilsonAPI>();
    std::shared_ptr<ICoreAPI<Model>> model_api = std::make_shared<ModelAPI>();
    std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api = std::make_shared<ScaleSetter>(ScaleType::MATCHING);
    std::shared_ptr<ICoreAPI<std::string>> marty_model_name = std::make_shared<MartyModelNameAPI>();
    std::shared_ptr<ICoreAPI<fs::path>> marty_model_path = std::make_shared<MartyModelPathAPI>();
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;
    std::shared_ptr<IMartyWilsonPathProxy> marty_paths = std::make_shared<MartyWilsonPathProxy>();
    std::shared_ptr<ICoreAPI<bool>> hard_coded_lo =
    std::make_shared<SMFromHypProxy>();
    if (use_marty->get()) {
        marty_proxy = std::make_shared<MartyWilsonProxy>(); 
    }
    WilsonGroupAdapterConfig adapters(wilson_proxy, iblock_c, use_marty, marty_model_name, marty_model_path, marty_proxy);
    this->current_group_adapters = std::make_shared<WilsonGroupAdapterConfig>(adapters);
    this->current_marty_paths = marty_paths;

    auto reg_ptr = make_registry();
    auto build_group_fn = [reg_ptr, adapters, marty_paths](WGroupId gid, Model mdl, bool useMarty, ContributionType ct, std::string group_name = "") -> std::shared_ptr<CoefficientGroup> {
        BuildContext ctx{
            .adapters = adapters,
            .model    = mdl,
            .backend  = useMarty ? Backend::Marty : Backend::Builtin,
            .contrib  = ct,
            .group_id = gid,
            .group_name = std::move(group_name),
            .marty_paths = marty_paths
        };
        CoefficientGroupBuilder b{*reg_ptr};
        return b.build(ctx);
    };

    std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
    const bool marty = use_marty->get();

    for (auto& g_id : config.groups) {
        ContributionType ct = (model == Model::SM) ? ContributionType::SM : ContributionType::BSM;
        auto grp = build_group_fn(g_id, model, marty, ct);
        groups.emplace(GroupMapper::str(g_id), std::move(grp));

    }
    if (use_marty->get()
        && config.order > QCDOrder::LO
        && !(hard_coded_lo && hard_coded_lo->get()))
    {
        // Marty pur => LO-only global
        config.order = QCDOrder::LO;
    }

    WilsonPortsConfig port_config{iblock_c, wilson_proxy, use_marty, has_wilson, model_api, scale_setter_api, hard_coded_lo, marty_paths};
    if (use_marty->get()) {
        port_config.matching_patches.push_back(make_hyperiso_c9_light_photon_sm_patch());
    }
    port_config.build_group = build_group_fn;

    this->cm = CoefficientManager::Builder(groups, config.matching_scale, config.hadronic_scale, OrderMapper::str(config.order), port_config, wilson_param_helpers);
}

void WilsonBuilder::add(WilsonBuildConfig config) {

    if (!this->cm) {
        LOG_ERROR("LogicError", "Cannot add Wilson groups before WilsonBuilder has been built.");
    }

    if (!this->current_group_adapters) {
        LOG_ERROR("LogicError", "WilsonBuilder has no cached group adapters for Wilson group registration.");
    }

    auto iblock_c = this->current_group_adapters->iblock_c;

    if (!wilson_param_helpers[Model::SM]) {
        wilson_param_helpers[Model::SM] = std::make_shared<WilsonParameterHelper>(iblock_c);
    }

    for (const auto& elem : config.groups) {
        wilson_param_helpers[Model::SM]->init(2, elem);
    }

    Model model = ModelAPI().get();

    if (model == Model::THDM) {
        if (!wilson_param_helpers[model]) {
            wilson_param_helpers[model] = std::make_shared<THDMParameterHelper>(iblock_c);
        }
        for (const auto& elem : config.groups) {
            wilson_param_helpers[model]->init(2, elem);
        }
    } else if (model == Model::SUSY) {
        if (!wilson_param_helpers[model]) {
            wilson_param_helpers[model] = std::make_shared<SUSYParameterHelper>(iblock_c);
        }
        for (const auto& elem : config.groups) {
            wilson_param_helpers[model]->init(2, elem);
        }
    }

    this->cm->set_matching_scale(config.matching_scale);
    this->cm->set_hadronic_scale(config.hadronic_scale);

    const bool marty = this->current_group_adapters->use_marty->get();
    auto hard_coded_lo = std::make_shared<SMFromHypProxy>();
    if (marty
        && config.order > QCDOrder::LO
        && !(hard_coded_lo && hard_coded_lo->get()))
    {
        config.order = QCDOrder::LO;
    }

    auto reg_ptr = make_registry();
    CoefficientGroupBuilder b{*reg_ptr};
    auto marty_paths = this->current_marty_paths
        ? this->current_marty_paths
        : std::make_shared<MartyWilsonPathProxy>();

    for (auto& g_id : config.groups) {
        ContributionType ct = (model == Model::SM) ? ContributionType::SM : ContributionType::BSM;

        BuildContext ctx{
            .adapters = *this->current_group_adapters,
            .model = model,
            .backend = marty ? Backend::Marty : Backend::Builtin,
            .contrib = ct,
            .group_id = g_id,
            .marty_paths = marty_paths
        };

        auto grp = b.build(ctx);

        std::string group_name = GroupMapper::str(g_id);
        LOG_VERBOSE("Initializing group", group_name);
        this->cm->registerCoefficientGroup(group_name, std::move(grp));
        LOG_VERBOSE("Initializing group at matching scale", group_name);
        this->cm->init_group_matching(group_name, OrderMapper::str(config.order));
        LOG_VERBOSE("Initializing group at hardronic scale", group_name);
        this->cm->init_group_hadronic_all_bases(group_name, OrderMapper::str(config.order));
        LOG_VERBOSE("Initialization done");
    }
}

void WilsonBuilder::add_custom_group(const CustomWilsonGroupConfig& config) {
    if (!this->cm) {
        LOG_ERROR("LogicError", "Cannot add a custom Wilson group before WilsonBuilder has been built.");
    }

    if (!this->current_group_adapters) {
        LOG_ERROR("LogicError", "WilsonBuilder has no cached group adapters for custom Wilson group registration.");
    }

    if (config.coefficients.empty()) {
        LOG_ERROR("ValueError", "Cannot add custom Wilson group", GroupMapper::str(config.group), "without coefficients.");
    }

    auto group = make_custom_wilson_group(config, *this->current_group_adapters);
    const std::string group_name = GroupMapper::str(config.group);

    this->cm->set_matching_scale(config.matching_scale);
    this->cm->set_hadronic_scale(config.hadronic_scale);
    this->cm->registerCoefficientGroup(group_name, std::move(group));
    this->cm->init_group_matching(group_name, OrderMapper::str(config.order));
    this->cm->init_group_hadronic_all_bases(group_name, OrderMapper::str(config.order));
}


void WilsonBuilder::add_matching_patch(const WilsonMatchingPatch& patch) {
    if (!this->cm) {
        LOG_ERROR("LogicError", "Cannot add a Wilson matching patch before WilsonBuilder has been built.");
    }
    this->cm->add_matching_patch(patch);
}

void WilsonBuilder::add_matching_patches(const std::vector<WilsonMatchingPatch>& patches) {
    if (!this->cm) {
        LOG_ERROR("LogicError", "Cannot add Wilson matching patches before WilsonBuilder has been built.");
    }
    this->cm->add_matching_patches(patches);
}

std::shared_ptr<WilsonProvider> WilsonBuilder::get_wilson_provider() {
    return std::make_shared<WilsonProvider>(shared_from_this());
}

std::shared_ptr<CoefficientManager> WilsonBuilder::get_coefficient_manager() {
    return cm;
}