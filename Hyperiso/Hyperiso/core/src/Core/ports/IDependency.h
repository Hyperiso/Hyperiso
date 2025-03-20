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
                                DepUpdateFunc recalculateFunc) 
                                = 0;

    virtual void remove_dependency(const std::string& name, ParameterType src) = 0;

    virtual void update_dependency(const std::string& name, ParameterType src) = 0;
};

#endif // __IDEPENDECY_H__
