#ifndef __WILSONFREEZER_H__
#define __WILSONFREEZER_H__

#include "IWilsonFreezer.h"
#include "Freezer.h"
#include "Include.h"
#include "ObsWilsonProxy.h"
#include "ObsWilsonBuilder.h"

class WilsonFreezer : public IWilsonFreezer<WGroup> {
public:
    WilsonFreezer(const std::shared_ptr<ObsWilsonBuilder> &wil_builder) {
        this->w_proxy = wil_builder->get_proxy();
    }

    void freeze(WGroup group) override;
    void unfreeze(WGroup group) override;

private:
    std::shared_ptr<ObsWilsonProxy> w_proxy;
};

#endif // __WILSONFREEZER_H__
