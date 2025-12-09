#ifndef STAT_PARAM_OPTIMIZER_PROXY_H
#define STAT_PARAM_OPTIMIZER_PROXY_H

#include "ParamOptimizerAdapter.h"
#include "IStatParamOptimizerProxy.h"

class StatParamOptimizerProxy : public IStatParamOptimizerProxy {
public:
    StatParamOptimizerProxy();
    
    void set_value(const BlockName& block, const LhaID& id, scalar_t v) override;
    void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) override;

    void remove(const BlockName& block, const LhaID& id) override;

    void commit(bool coalesce = true) override;

    void clear() override;

private:
    ParamOptimizerAdapter poa;

};

#endif