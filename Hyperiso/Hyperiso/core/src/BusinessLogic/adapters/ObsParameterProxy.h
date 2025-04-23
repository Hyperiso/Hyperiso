#ifndef MODEL_PARAM_ADAPTER_H
#define MODEL_PARAM_ADAPTER_H

#include "IObsParameterAdapter.h"
#include "ParameterProvider.h"
#include "Include.h"

class ObsParameterProxy : public IObsParameterProxy<std::string, LhaID> {
public:
ObsParameterProxy(ParameterType type);

    scalar_t operator()(const ParamId& pid, ParameterProvider::DataType d_type=ParameterProvider::DataType::VALUE);
    scalar_t operator()(const std::string& block, const LhaID& id) const;
    
private:
    ParameterProvider pp;
    ParameterProvider pp_with_type;
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON};
};

#endif 