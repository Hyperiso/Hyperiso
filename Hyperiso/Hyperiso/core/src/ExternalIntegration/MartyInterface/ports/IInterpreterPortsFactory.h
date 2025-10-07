#ifndef IINTERPRETER_PORTS_FACTORY_H
#define IINTERPRETER_PORTS_FACTORY_H

#include <memory>
#include <string>
#include "IParameterResolver.h"

class IInterpreterPortsFactory {
public:
    virtual ~IInterpreterPortsFactory() = default;

    virtual std::unique_ptr<IParameterResolver>
    makeResolver(const std::string& modelName,
                 const std::string& modelJsonPath,
                 const std::string& smJsonPath) const = 0;
};

#endif 
