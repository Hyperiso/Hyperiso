#include "ObservableInterface.h"

ObservableInterface::ObservableInterface() {
    WilsonBuildConfig config;
    auto builder = std::make_shared<ObsWilsonBuilder>(config);
    manager = std::make_shared<ObsManager>(builder);
}