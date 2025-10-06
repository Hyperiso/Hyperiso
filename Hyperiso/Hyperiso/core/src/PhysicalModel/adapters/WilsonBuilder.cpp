#include "WilsonBuilder.h"

WilsonBuilder::WilsonBuilder(WilsonBuildConfig config) {
    this->build(config);
}

WilsonBuilder::WilsonBuilder(std::shared_ptr<CoefficientManager> manager) : cm(manager) {}

void WilsonBuilder::build(WilsonBuildConfig config) {
    // WilsonParameterHelper().init(2);
    std::map<Model, std::shared_ptr<IWilsonParameterHelper>> wilson_param_helpers {};

    std::shared_ptr<IBlockComposer> iblock_c = std::make_shared<WilsonParamComposer>();
    
    wilson_param_helpers[Model::SM] = std::make_shared<WilsonParameterHelper>(iblock_c);
    wilson_param_helpers[Model::SM]->init(2);
    std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
    Model model = ModelAPI().get();
    if (model == Model::THDM) {
        wilson_param_helpers[model] = std::make_shared<thdm_parameters>(iblock_c);
        wilson_param_helpers[model]->init(2);
    } else if (model == Model::SUSY) {
        wilson_param_helpers[model] = std::make_shared<susy_parameters>(iblock_c);
        wilson_param_helpers[model]->init(2);
    }
    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy = std::make_shared<ParameterProxy>(ParameterType::WILSON);
    std::shared_ptr<IParameterProxy<std::string, LhaID>> sm_proxy = std::make_shared<ParameterProxy>(ParameterType::SM);
    std::shared_ptr<ICoreAPI<bool>> use_marty = std::make_shared<UseMarty>();
    std::shared_ptr<ICoreAPI<Model>> model_api = std::make_shared<ModelAPI>();
    std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api = std::make_shared<ScaleSetter>(ScaleType::MATCHING);
    std::shared_ptr<ICoreAPI<std::string>> marty_model_name = std::make_shared<MartyModelNameAPI>();
    std::shared_ptr<ICoreAPI<fs::path>> marty_model_path = std::make_shared<MartyModelPathAPI>();

    WilsonGroupAdapterConfig adapters(wilson_proxy, iblock_c, use_marty, marty_model_name, marty_model_path);

    for (auto& g_id : config.groups) {
        std::string gn_str = GroupMapper::str(g_id);
        groups.emplace(gn_str, WilsonGroupFactory::create_coefficient_group(g_id, model, adapters));
    }
    std::cout << "FUC3K " << std::endl; 
    if (UseMarty().get() && config.order > QCDOrder::LO) {
        // TODO : Use QCDOrder > LO as a flag for SM calculation. 
        LOG_WARN("Using MARTY defaults all calculations to LO.");
        config.order = QCDOrder::LO;
    }

    PortsConfig port_config{iblock_c, wilson_proxy, use_marty, model_api, scale_setter_api};

    this->cm = CoefficientManager::Builder(ModelMapper::str(model), groups, config.matching_scale, config.hadronic_scale, OrderMapper::str(config.order), port_config, wilson_param_helpers);
}

void WilsonBuilder::add(WilsonBuildConfig config) {
    Model model = ModelAPI().get();
    this->cm->set_matching_scale(config.matching_scale);
    this->cm->set_hadronic_scale(config.hadronic_scale);

    std::shared_ptr<IBlockComposer> iblock_c = std::make_shared<WilsonParamComposer>();
    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy = std::make_shared<ParameterProxy>(ParameterType::WILSON);
    std::shared_ptr<IParameterProxy<std::string, LhaID>> sm_proxy = std::make_shared<ParameterProxy>(ParameterType::SM);
    std::shared_ptr<ICoreAPI<bool>> use_marty = std::make_shared<UseMarty>();
    std::shared_ptr<ICoreAPI<Model>> model_api = std::make_shared<ModelAPI>();
    std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api = std::make_shared<ScaleSetter>(ScaleType::MATCHING);
    std::shared_ptr<ICoreAPI<std::string>> marty_model_name = std::make_shared<MartyModelNameAPI>();
    std::shared_ptr<ICoreAPI<fs::path>> marty_model_path = std::make_shared<MartyModelPathAPI>();

    WilsonGroupAdapterConfig adapters(wilson_proxy, iblock_c, use_marty, marty_model_name, marty_model_path);

    for (auto& g_id : config.groups) {
        std::string group_name = GroupMapper::str(g_id);
        LOG_INFO("Initializing group", group_name);
        this->cm->registerCoefficientGroup(group_name, WilsonGroupFactory::create_coefficient_group(g_id, model, adapters));
        this->cm->init_group_matching(group_name, OrderMapper::str(config.order));
        this->cm->init_group_hadronic_all_bases(group_name, OrderMapper::str(config.order));
    }
}

std::shared_ptr<WilsonProvider> WilsonBuilder::get_wilson_provider() {
    return std::make_shared<WilsonProvider>(shared_from_this());
}

std::shared_ptr<CoefficientManager> WilsonBuilder::get_coefficient_manager() {
    return cm;
}