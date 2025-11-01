#ifndef USE_MARTY_H
#define USE_MARTY_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

class UseMarty : public ICoreAPI<bool> {
public:
    inline bool get() override {return HyperisoMaster().get_model() == Model::MARTY ? 1 : 0;}
};

#endif