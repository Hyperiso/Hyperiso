#ifndef __COMPOSITEPARAMCREATOR_H__
#define __COMPOSITEPARAMCREATOR_H__

#include "IDependency.h"
#include "DependentBlockManager.h"

class CompositeParamAdapter : public IDependency {
private:
    /* data */
public:
    void add_block_dependency(
        const std::string& name,
        const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
        ParameterType dest,
        DepUpdateFunc recalculateFunc
    ) override;

    void add_param_dependency(
        const ParamId& pid,
        const std::unordered_set<ParamId>& source_pids,
        DepParamUpdateFunc recalculateFunc
    ) override;

    void remove_dependency(const std::string& name, ParameterType src) override;

    void update_dependency(const std::string& name, ParameterType src) override;

};

#endif // __COMPOSITEPARAMCREATOR_H__
