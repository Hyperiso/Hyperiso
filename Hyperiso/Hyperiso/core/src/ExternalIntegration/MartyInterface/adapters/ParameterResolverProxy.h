#ifndef PARAMETER_RESOLVER_PROXY_H
#define PARAMETER_RESOLVER_PROXY_H

#include <unordered_map>
#include <string>
#include <memory>
#include "IParameterResolver.h"
#include "IMappingDatabasePort.h"
#include "Logger.h"
#include <iomanip>

class ParameterResolverProxy final : public IParameterResolver {
public:
    ParameterResolverProxy(std::shared_ptr<IMappingDatabasePort> sm,
                           std::shared_ptr<IMappingDatabasePort> model);

    std::unordered_map<std::string, ResolvedParam>
    resolve(const std::vector<Extractor::Parameter>& params,
            bool modelIsSM) const override;

    std::unique_ptr<IParameterResolver> clone() const override;

private:
    std::shared_ptr<IMappingDatabasePort> sm_;
    std::shared_ptr<IMappingDatabasePort> model_;
};

#endif // PARAMETER_RESOLVER_PROXY_H
