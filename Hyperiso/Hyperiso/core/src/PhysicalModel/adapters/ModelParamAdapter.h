#ifndef MODEL_PARAM_ADAPTER_H
#define MODEL_PARAM_ADAPTER_H

#include "IParamAdapter.h"
#include "ParameterProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

class ParameterProxy : public IParameterProxy<std::string, LhaID> {
public:
    ParameterProxy(ParameterType type);
    double operator()(const std::string& block, const LhaID& id) override;
    
private:
    ParameterProvider pp;
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::THDM, ParameterType::SUSY, ParameterType::WILSON};
};

#endif 