#include "ParameterResolverProxy.h"

ParameterResolverProxy::ParameterResolverProxy(std::shared_ptr<IMappingDatabasePort> sm,
                        std::shared_ptr<IMappingDatabasePort> model)
    : sm_(std::move(sm)), model_(std::move(model)) {}

std::unordered_map<std::string, ResolvedParam>
ParameterResolverProxy::resolve(const std::vector<Extractor::Parameter>& params,
        bool modelIsSM) const
{
    std::unordered_map<std::string, ResolvedParam> out;
    auto smMap    = sm_->getParams();
    auto modelMap = model_->getParams();
    
    for (const auto& p : params) {
        const auto& name = p.name;
        auto it = smMap.find(name);
        if (it != smMap.end()) {
            
            out[name] = { it->second.block, it->second.code, p.complex, false };
            continue;
        }
        it = modelMap.find(name);
        if (it != modelMap.end()) {
            if (modelIsSM) {
                LOG_ERROR("LogicError",
                    "Accès à un paramètre BSM pendant un calcul SM. Vérifie les fichiers de mapping MARTY.");
            }
            out[name] = { it->second.block, it->second.code, p.complex, true };
        } else {
            LOG_ERROR("Mapping", std::string("Paramètre introuvable: ") + name);
        }
    }
    return out;
}

std::unique_ptr<IParameterResolver> ParameterResolverProxy::clone() const {
    return std::make_unique<ParameterResolverProxy>(sm_, model_);
}