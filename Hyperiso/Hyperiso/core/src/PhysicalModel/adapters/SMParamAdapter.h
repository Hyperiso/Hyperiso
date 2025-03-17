#ifndef SM_PARAM_ADAPTER_H
#define SM_PARAM_ADAPTER_H

#include <string>
#include "IParamAdapter.h"
#include "ParameterProvider.h"
#include "General.h"

class SMParamAdapter : public IParamAdapter<std::string, int> {
public:
        double operator()(std::string block, int code) override;

private:
    ParameterProvider pp{ParameterType::SM};
};

#endif 