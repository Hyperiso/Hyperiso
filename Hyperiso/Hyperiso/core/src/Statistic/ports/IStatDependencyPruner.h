#ifndef __ISTATDEPENDENCYPRUNER_H__
#define __ISTATDEPENDENCYPRUNER_H__

#include "Include.h"

class IStatDependencyPruner {
public:
    virtual ~IStatDependencyPruner() = default;

    virtual void reattach_block(ParameterType, const std::string&) = 0;
    virtual void detach_block(ParameterType, const std::string&) = 0;

    virtual void reattach_parameter(ParameterType, const std::string&, const LhaID&) = 0;
    virtual void detach_parameter(ParameterType, const std::string&, const LhaID&) = 0;
};

#endif // __ISTATDEPENDENCYPRUNER_H__
