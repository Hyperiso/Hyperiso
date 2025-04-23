#ifndef USE_MARTY_H
#define USE_MARTY_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "IObsCoreAPI.h"

class UseMarty : public IObsCoreAPI<bool> {
public:
    inline bool get() override {return HyperisoMaster().check_flag(ExternalFlag::USE_MARTY);}
};

#endif