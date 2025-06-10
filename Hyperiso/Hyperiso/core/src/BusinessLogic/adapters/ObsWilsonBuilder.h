#ifndef OBS_WILSON_BUILDER_H
#define OBS_WILSON_BUILDER_H

#include "IObsWilsonBuilder.h"
#include "IWilsonBuilder.h"
#include "ObsWilsonProxy.h"
#include "WilsonProvider.h"
#include "Configs.h"
#include "Include.h"

class ObsWilsonBuilder : public IObsWilsonBuilder<ObsWilsonProxy, WGroup> {
public:
    ObsWilsonBuilder(std::shared_ptr<WilsonBuilder> wil_builder) : wil_builder(std::move(wil_builder)) {}

    void build(std::shared_ptr<AbstractConfig> config) override;
    std::shared_ptr<ObsWilsonProxy> get_proxy() override;

private:
    std::shared_ptr<WilsonBuilder> wil_builder;
};

#endif