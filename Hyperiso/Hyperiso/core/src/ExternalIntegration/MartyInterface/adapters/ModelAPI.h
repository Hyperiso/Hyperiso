#ifndef MODEL_API_H
#define MODEL_API_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

class ModelAPI : public ICoreAPI<Model> {
public:
    inline Model get() override {return HyperisoMaster().get_model();}
};

#endif