#ifndef __COMPOSITEPARAMCREATOR_H__
#define __COMPOSITEPARAMCREATOR_H__

#include "IDependency.h"
#include "DependentBlockManager.h"

class CompositeParamCreator : public IDependency {
private:
    /* data */
public:
    void add_dependency(
        const std::string& name,
        const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
        ParameterType dest,
        std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> recalculateFunc
    ) override;

    void remove_dependency(const std::string& name, ParameterType src) override;

};

#endif // __COMPOSITEPARAMCREATOR_H__
