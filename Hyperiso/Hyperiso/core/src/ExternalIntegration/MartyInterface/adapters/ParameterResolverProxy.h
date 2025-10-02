// ParameterResolverProxy.h
#ifndef PARAMETER_RESOLVER_PROXY_H
#define PARAMETER_RESOLVER_PROXY_H

#include <unordered_map>
#include <string>
#include <memory>
#include "IParameterResolver.h"
#include "IMappingDatabasePort.h"
#include "Logger.h"

class ParameterResolverProxy final : public IParameterResolver {
public:
    ParameterResolverProxy(std::shared_ptr<IMappingDatabasePort> sm,
                           std::shared_ptr<IMappingDatabasePort> model)
        : sm_(std::move(sm)), model_(std::move(model)) {}

    std::unordered_map<std::string, ResolvedParam>
    resolve(const std::vector<Extractor::Parameter>& params,
            bool modelIsSM) const override
    {
        std::unordered_map<std::string, ResolvedParam> out;
        const auto& smMap    = sm_->getParams();
        const auto& modelMap = model_->getParams();

        for (const auto& p : params) {
            const auto& name = p.name;
            auto it = smMap.find(name);
            if (it != smMap.end()) {
                out[name] = { it->second.block, it->second.pdgCode, p.complex, false };
                continue;
            }
            it = modelMap.find(name);
            if (it != modelMap.end()) {
                if (modelIsSM) {
                    LOG_ERROR("LogicError",
                        "Accès à un paramètre BSM pendant un calcul SM. Vérifie les fichiers de mapping MARTY.");
                }
                out[name] = { it->second.block, it->second.pdgCode, p.complex, true };
            } else {
                LOG_ERROR("Mapping", std::string("Paramètre introuvable: ") + name);
            }
        }
        return out;
    }

    std::unique_ptr<IParameterResolver> clone() const override {
        return std::make_unique<ParameterResolverProxy>(sm_, model_);
    }

private:
    std::shared_ptr<IMappingDatabasePort> sm_;
    std::shared_ptr<IMappingDatabasePort> model_;
};

#endif // PARAMETER_RESOLVER_PROXY_H
