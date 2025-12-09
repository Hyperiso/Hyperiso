#ifndef PARAM_OPTIMIZER_ADAPTER_H
#define PARAM_OPTIMIZER_ADAPTER_H

#include "ParamOptimizer.h"
#include "IParamOptimizer.h"

class ParamOptimizerAdapter : public IParamOptimizer {
public:
    ParamOptimizerAdapter(std::vector<ParameterType> scopes);
    
    void set_value(const BlockName& block, const LhaID& id, scalar_t v) override;
    void set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) override;

    void remove(const BlockName& block, const LhaID& id) override;

    void commit(bool coalesce = true) override;

    void clear() override;

private:
    std::shared_ptr<ParamOptimizer> po;

};

#endif