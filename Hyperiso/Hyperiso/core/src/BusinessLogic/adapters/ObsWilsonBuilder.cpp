#include "ObsWilsonBuilder.h"

void ObsWilsonBuilder::build(std::shared_ptr<AbstractConfig> config) {
    auto wil_config = *std::dynamic_pointer_cast<WilsonBuildConfig>(config);

    if (built) {
        wil_builder->add(wil_config);
    } else {
        wil_builder->build(wil_config);
        built = true;
    }
}

std::shared_ptr<ObsWilsonProxy> ObsWilsonBuilder::get_proxy() {
    return std::make_shared<ObsWilsonProxy>(wil_builder->get_wilson_provider());
}

void ObsWilsonBuilder::switch_basis(WGroup group) {
    this->wil_builder->switch_basis(group);
}
