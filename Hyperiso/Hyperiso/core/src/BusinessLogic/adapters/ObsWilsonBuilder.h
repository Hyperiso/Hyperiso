#ifndef OBS_WILSON_BUILDER_H
#define OBS_WILSON_BUILDER_H

#include "IObsWilsonBuilder.h"
#include "IWilsonBuilder.h"
#include "ObsWilsonProxy.h"
#include "Configs.h"
#include "Include.h"

class ObsWilsonBuilder : public IObsWilsonBuilder<ObsWilsonProxy, WGroup> {
public:
    ObsWilsonBuilder() : wil_builder(std::make_shared<WilsonBuilder>()) {}
    ObsWilsonBuilder(std::shared_ptr<IWilsonBuilder<WilsonBuildConfig, WGroup, WilsonProvider>> wil_builder) : wil_builder(wil_builder) {}

    void build(std::shared_ptr<AbstractConfig> config) override;
    void switch_basis(WGroup group) override;
    std::shared_ptr<ObsWilsonProxy> get_proxy() override;

private:
    std::shared_ptr<IWilsonBuilder<WilsonBuildConfig, WGroup, WilsonProvider>> wil_builder;
    bool built {false};
};

#endif