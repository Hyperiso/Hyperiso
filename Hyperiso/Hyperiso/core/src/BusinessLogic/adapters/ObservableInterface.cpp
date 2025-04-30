#include "ObservableInterface.h"

ObservableInterface::ObservableInterface() {
    auto builder = std::make_shared<ObsWilsonBuilder>();
    manager = std::make_shared<ObsManager>(builder);
}