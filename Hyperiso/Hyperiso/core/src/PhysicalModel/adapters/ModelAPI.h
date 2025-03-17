#ifndef HAS_WILSON_API_H
#define HAS_WILSON_API_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

class HAS_WILSON_API : public ICoreAPI<Model> {
public:
    inline Model get() override {HyperisoMaster().get_model();}
};

#endif