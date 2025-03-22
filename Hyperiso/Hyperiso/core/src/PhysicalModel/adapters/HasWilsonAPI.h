#ifndef HAS_WILSON_API_H
#define HAS_WILSON_API_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

class HAS_WILSON_API : public ICoreAPI<bool> {
public:
    inline bool get() override {return HyperisoMaster().check_flag(ExternalFlag::HAS_WILSON_INPUT);}
};

#endif