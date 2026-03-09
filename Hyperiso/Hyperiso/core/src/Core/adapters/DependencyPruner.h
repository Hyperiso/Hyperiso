#ifndef __DEPENDENCYPRUNER_H__
#define __DEPENDENCYPRUNER_H__

#include "IDependencyPruner.h"
#include "BlockName.h"
#include "Include.h"
#include "Parameters.h"

class DependencyPruner : public IDependencyPruner<ParameterType, BlockName, LhaID> {
public:
    void reattach_block(ParameterType tp, const BlockName& block_name) override;
    void detach_block(ParameterType tp, const BlockName& block_name) override;

    void reattach_parameter(ParameterType tp, const BlockName& block_name, const LhaID& id) override;
    void detach_parameter(ParameterType tp, const BlockName& block_name, const LhaID& id) override;
};

#endif // __DEPENDENCYPRUNER_H__
