#ifndef ISTAT_PARAMETER_PROXY_H
#define ISTAT_PARAMETER_PROXY_H


#include "CorrelationProvider.h"
#include "Include.h"

struct IStatParameterProxy {
    virtual ~IStatParameterProxy() = default;

    virtual std::shared_ptr<Parameter> get_param(const ParamId&) const = 0;
    virtual std::shared_ptr<Parameter> get_param(const std::string& block, const LhaID& id) const = 0;
    virtual scalar_t operator()(const ParamId&, DataType d_type=DataType::VALUE) const = 0;
    virtual double operator()(const ObservableId&, DataType d_type=DataType::VALUE) const = 0;
    virtual scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const = 0 ;
    virtual std::shared_ptr<Parameter> get_obs_param(const ObservableId&) const = 0;
};

#endif