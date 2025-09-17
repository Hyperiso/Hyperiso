#include "ObservableInterface.h"

ObservableInterface::ObservableInterface() {
    WilsonBuildConfig config;
    std::shared_ptr<WilsonBuilder> builder_ptr = std::make_shared<WilsonBuilder>(config);
    auto builder = std::make_shared<ObsWilsonBuilder>(builder_ptr);
    manager = std::make_shared<ObsManager>(builder);
}

void ObservableInterface::add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr) {
    manager->add_custom_decay(id, ptr);
}