#ifndef __WILSONFREEZER_H__
#define __WILSONFREEZER_H__

#include "IWilsonFreezer.h"
#include "Freezer.h"
#include "Include.h"
#include "ObsWilsonProxy.h"
#include "IObsWilsonBuilder.h"

class WilsonFreezer : public IWilsonFreezer<WGroup> {
public:
    WilsonFreezer(std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> wil_builder) : w_proxy(wil_builder->get_proxy()) {}

    void freeze(WGroup group) override;
    void unfreeze(WGroup group) override;

private:
    std::shared_ptr<IObsWilsonProxy<ObsWilsonBuilder>> w_proxy;
};

#endif // __WILSONFREEZER_H__
