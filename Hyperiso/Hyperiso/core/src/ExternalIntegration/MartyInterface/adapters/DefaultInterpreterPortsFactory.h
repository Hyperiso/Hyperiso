#ifndef DEFAULT_INTERPRETER_PORTS_FACTORY_H
#define DEFAULT_INTERPRETER_PORTS_FACTORY_H

#include "IInterpreterPortsFactory.h"
#include <memory>

class IMappingAdapterFactory;
class IParamMappingSource;

class DefaultInterpreterPortsFactory : public IInterpreterPortsFactory {
public:
    explicit DefaultInterpreterPortsFactory(
        std::shared_ptr<IMappingAdapterFactory> adapterFactory = nullptr,
        std::shared_ptr<IParamMappingSource>    loader         = nullptr
    );

    std::unique_ptr<IParameterResolver>
    makeResolver(const std::string& modelName,
                 const std::string& modelJsonPath,
                 const std::string& smJsonPath) const override;

private:
    std::shared_ptr<IMappingAdapterFactory> adapterFactory_;
    std::shared_ptr<IParamMappingSource>    loader_;
};

#endif
