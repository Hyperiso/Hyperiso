#ifndef PARAMETER_PROXY_H
#define PARAMETER_PROXY_H

#include "IParameterProxy.h"
#include "ParameterProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

class ParameterProxy : public IParameterProxy<std::string, LhaID> {
public:
    ParameterProxy(ParameterType type);
    scalar_t operator()(const std::string& block, const LhaID& id) const override;
    
    bool exist(const std::string& block, const LhaID& id) const override;

    double get_scale(const std::string& block) const override;
    
private:
    ParameterProvider pp;
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON};
};

#endif 