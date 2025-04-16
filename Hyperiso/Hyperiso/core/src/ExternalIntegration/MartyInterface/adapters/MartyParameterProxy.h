#ifndef __MARTYPARAMETERPROXY_H__
#define __MARTYPARAMETERPROXY_H__

#include "IMartyParameterProxy.h"
#include "ParameterProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

class MartyParameterProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    MartyParameterProxy(ParameterType type);
    scalar_t operator()(const std::string& block, const LhaID& id) const override;
    
private:
    ParameterProvider pp;
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM};
};

#endif // __MARTYPARAMETERPROXY_H__
