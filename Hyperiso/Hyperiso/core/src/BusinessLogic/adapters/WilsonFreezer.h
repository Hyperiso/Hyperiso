#ifndef __WILSONFREEZER_H__
#define __WILSONFREEZER_H__

#include "IWilsonFreezer.h"
#include "Freezer.h"
#include "Include.h"

class WilsonFreezer : public IWilsonFreezer<WGroup> {
public:
    void freeze(WGroup group) override;
    void unfreeze(WGroup group) override;
};

#endif // __WILSONFREEZER_H__
