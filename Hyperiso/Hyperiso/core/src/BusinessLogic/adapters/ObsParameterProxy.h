#ifndef OBS_PARAMETER_PROXY_H
#define OBS_PARAMETER_PROXY_H

#include "IObsParameterProxy.h"
#include "ParameterProvider.h"
#include "Include.h"

class ObsParameterProxy : public IObsParameterProxy<std::string, LhaID> {
public:
    ObsParameterProxy(ParameterType type = ParameterType::SM);

    scalar_t operator()(const ParamId& pid, ParameterProvider::DataType d_type=ParameterProvider::DataType::VALUE);
    scalar_t operator()(const std::string& block, const LhaID& id, ParameterProvider::DataType d_type=ParameterProvider::DataType::VALUE) const;
    
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const;
private:
    ParameterProvider pp;
    ParameterProvider pp_with_type;
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON, ParameterType::FLAVOR};
};

#endif 