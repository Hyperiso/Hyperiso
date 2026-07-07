#include "ObsWilsonBuilder.h"

void ObsWilsonBuilder::build(std::shared_ptr<AbstractConfig> config) {
    auto wil_config = *std::dynamic_pointer_cast<WilsonBuildConfig>(config);

    if (wil_builder->get_coefficient_manager()) {
        wil_builder->add(wil_config);
    } else {
        wil_builder->build(wil_config);
    }
}


void ObsWilsonBuilder::add_custom_group(const CustomWilsonGroupConfig& config) {
    if (!wil_builder->get_coefficient_manager()) {
        WilsonBuildConfig base;
        base.matching_scale = config.matching_scale;
        base.hadronic_scale = config.hadronic_scale;
        base.order = config.order;
        wil_builder->build(base);
    }

    wil_builder->add_custom_group(config);
}


void ObsWilsonBuilder::add_matching_patch(const WilsonMatchingPatch& patch) {
    if (!wil_builder->get_coefficient_manager()) {
        WilsonBuildConfig base;
        base.order = patch.order == QCDOrder::NONE ? QCDOrder::LO : patch.order;
        wil_builder->build(base);
    }

    wil_builder->add_matching_patch(patch);
}

std::shared_ptr<IObsWilsonProxy> ObsWilsonBuilder::get_proxy() {
    return std::make_shared<ObsWilsonProxy>(wil_builder->get_wilson_provider());
}
