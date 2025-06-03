#include "WilsonBuilder.h"

//DUPLICATE the creation of coefficientmanager -> lead to ~CoefficientManager, destroy everything
// WilsonBuilder::WilsonBuilder() {
//     this->cm = std::make_shared<CoefficientManager>();
// }

WilsonBuilder::WilsonBuilder(WilsonBuildConfig config) {
    this->build(config);
}

WilsonBuilder::WilsonBuilder(std::shared_ptr<CoefficientManager> manager) : cm(manager) {}

void WilsonBuilder::build(WilsonBuildConfig config) {
    WilsonParameterHelper().init(2);
    std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
    Model model = ModelAPI().get();
    for (auto& g_id : config.groups) {
        std::string gn_str = GroupMapper::str(g_id);
        LOG_INFO("Creating group", gn_str);
        groups.emplace(gn_str, WilsonGroupFactory::create_coefficient_group(g_id, model));
    }

    if (UseMarty().get() && config.order > QCDOrder::LO) {
        LOG_WARN("Using MARTY defaults all calculations to LO.");
        config.order = QCDOrder::LO;
    }

    this->cm = CoefficientManager::Builder(ModelMapper::str(model), groups, config.matching_scale, config.hadronic_scale, OrderMapper::str(config.order));
}

void WilsonBuilder::add(WilsonBuildConfig config) {
    auto init_group = [&](WGroup g_id, Model model) {
        std::string group_name = GroupMapper::str(g_id);
        this->cm->registerCoefficientGroup(group_name, WilsonGroupFactory::create_coefficient_group(g_id, model));
        this->cm->init_group_matching(group_name, OrderMapper::str(config.order));
        this->cm->init_group_hadronic_all_bases(group_name, OrderMapper::str(config.order));
    };

    Model model = ModelAPI().get();

    for (auto& g_id : config.groups) {
        init_group(g_id, Model::SM);
        if (model == Model::THDM || model == Model::SUSY) {
            init_group(g_id, model);
        }
    }
}

std::shared_ptr<WilsonProvider> WilsonBuilder::get_wilson_provider() {
    return std::make_shared<WilsonProvider>(this->cm);
}
