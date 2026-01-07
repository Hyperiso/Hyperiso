#include "ObsWilsonBuilder.h"

void ObsWilsonBuilder::build(std::shared_ptr<AbstractConfig> config) {
    auto wil_config = *std::dynamic_pointer_cast<WilsonBuildConfig>(config);

    if (wil_builder->get_coefficient_manager()) {
        wil_builder->add(wil_config);
    } else {
        wil_builder->build(wil_config);
    }
}

std::shared_ptr<IObsWilsonProxy> ObsWilsonBuilder::get_proxy() {
    return std::make_shared<ObsWilsonProxy>(wil_builder->get_wilson_provider());
}
