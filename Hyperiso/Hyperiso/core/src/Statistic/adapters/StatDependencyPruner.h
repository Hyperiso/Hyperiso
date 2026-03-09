#ifndef __STATDEPENDENCYPRUNER_H__
#define __STATDEPENDENCYPRUNER_H__

#include "IStatDependencyPruner.h"
#include "DependencyPruner.h"

class StatDependencyPruner : public IStatDependencyPruner {
public:
    void reattach_block(ParameterType tp, const std::string& block_name) override;
    void detach_block(ParameterType tp, const std::string& block_name) override;

    void reattach_parameter(ParameterType tp, const std::string& block_name, const LhaID& id) override;
    void detach_parameter(ParameterType tp, const std::string& block_name, const LhaID& id) override;

private:
    DependencyPruner dp;
};

#endif // __STATDEPENDENCYPRUNER_H__
