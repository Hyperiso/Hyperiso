#ifndef __WILSONFREEZER_H__
#define __WILSONFREEZER_H__

#include "IWilsonFreezer.h"
#include "Freezer.h"
#include "Include.h"
#include "ObsWilsonProxy.h"
#include "ObsWilsonBuilder.h"

class WilsonFreezer : public IWilsonFreezer<WGroupId> {
public:
    WilsonFreezer(const std::shared_ptr<ObsWilsonBuilder> &wil_builder) {
        this->w_proxy = wil_builder->get_proxy();
    }

    void freeze(WGroupId group) override;
    void unfreeze(WGroupId group) override;

private:
    std::shared_ptr<ObsWilsonProxy> w_proxy;
};

#endif // __WILSONFREEZER_H__
