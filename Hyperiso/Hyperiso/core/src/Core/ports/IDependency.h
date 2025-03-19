#ifndef __IDEPENDECY_H__
#define __IDEPENDECY_H__

#include "Include.h"
#include "Block.h"

class IDependency {
public:
    virtual ~IDependency() = default;

    virtual void add_dependency(const std::string& name,
                                const std::unordered_map<ParameterType, std::vector<std::string>>& source_names,
                                ParameterType dest,
                                std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> recalculateFunc) 
                                = 0;

    virtual void remove_dependency(const std::string& name, ParameterType src) = 0;
};

#endif // __IDEPENDECY_H__
