#ifndef OBS_PARAMETER_PROXY_H
#define OBS_PARAMETER_PROXY_H

#include "IStatParameterProxy.h"
#include "ParameterProvider.h"
#include "Include.h"

class StatParameterProxy : public IStatParameterProxy {
public:
    StatParameterProxy(ParameterType type = ParameterType::SM);

    std::shared_ptr<Parameter> get_param(const ParamId&) const;
    std::shared_ptr<Parameter> get_param(const std::string& block, const LhaID& id) const;
    scalar_t operator()(const ParamId&, DataType d_type=DataType::VALUE) const;
    double operator()(const ObservableId&, DataType d_type=DataType::VALUE) const;
    scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const;
    std::shared_ptr<Parameter> get_obs_param(const ObservableId&) const;
private:
    ParameterProvider pp;
    ParameterProvider pp_with_type;
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON, ParameterType::FLAVOR, ParameterType::DECAY, ParameterType::OBSERVABLE};
};

#endif 