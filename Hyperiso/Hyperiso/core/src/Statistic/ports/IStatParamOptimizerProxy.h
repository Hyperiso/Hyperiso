#ifndef ISTAT_PARAM_OPTIMIZER_PROXY_H
#define ISTAT_PARAM_OPTIMIZER_PROXY_H

#include "Include.h"
#include "Parameter.h"

class IStatParamOptimizerProxy {
public:
    virtual ~IStatParamOptimizerProxy() = default;
    
    virtual void set_value(const BlockName& block, const LhaID& id, scalar_t v) = 0;
    virtual void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p)  = 0;

    virtual void remove(const BlockName& block, const LhaID& id)  = 0;

    virtual void commit(bool coalesce = true)  = 0;

    virtual void clear()  = 0;

};

#endif