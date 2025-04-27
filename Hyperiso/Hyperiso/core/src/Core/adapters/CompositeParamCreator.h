#ifndef COMPOSITEPARAMCREATOR_H
#define COMPOSITEPARAMCREATOR_H

#include "IDependency.h"
#include "DependentBlockManager.h"

/**
 * @class CompositeParamAdapter
 * @ingroup DependencyManagementModule
 * @brief Concrete implementation of IDependency using DependentBlockManager.
 *
 * This class adapts high-level operations to the underlying DependentBlockManager.
 */
class CompositeParamAdapter : public IDependency {
private:
    /* data */
public:

    /**
     * @brief Adds a dependent block using DependentBlockManager.
     */
    void add_block_dependency(
        const std::string& name,
        const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
        ParameterType dest,
        DepUpdateFunc recalculateFunc
    ) override;

    /**
     * @brief Adds a dependent parameter using DependentBlockManager.
     */
    void add_param_dependency(
        const ParamId& pid,
        const std::unordered_set<ParamId>& source_pids,
        DepParamUpdateFunc recalculateFunc
    ) override;

    /**
     * @brief Removes a dependent block using DependentBlockManager.
     */
    void remove_dependency(const std::string& name, ParameterType src) override;

    /**
     * @brief Updates a dependent block using DependentBlockManager.
     */
    void update_dependency(const std::string& name, ParameterType src) override;

};

#endif // COMPOSITEPARAMCREATOR_H
