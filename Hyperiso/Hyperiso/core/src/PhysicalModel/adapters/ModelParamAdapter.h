#ifndef MODEL_PARAM_ADAPTER_H
#define MODEL_PARAM_ADAPTER_H

#include <string>
#include "IParamAdapter.h"
#include "ParameterProvider.h"
#include "HyperisoMaster.h"
#include "General.h"

class ModelParamAdapter : public IParamAdapter<std::string, int> {
public:
        double operator()(std::string block, int code) override;
private:
    ParameterProvider pp{ParameterTypeMapper::enum_elt(ModelMapper::str(HyperisoMaster().get_model()))};
};

#endif 