#ifndef OBS_USE_MARTY_H
#define OBS_USE_MARTY_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "IObsCoreAPI.h"

class ObsUseMarty : public IObsCoreAPI<bool> {
public:
    inline bool get() override {return HyperisoMaster().get_model() == Model::MARTY ? 1 : 0;}
};

#endif