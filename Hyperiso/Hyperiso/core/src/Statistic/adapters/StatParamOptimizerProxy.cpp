#include "StatParamOptimizerProxy.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

#include "MemoryManager.h"


StatParamOptimizerProxy::StatParamOptimizerProxy()
    : poa_bsm(nullptr),
      poa_standard({
          ParameterType::SM,
          ParameterType::FLAVOR,
          ParameterType::DECAY,
          ParameterType::WILSON
      }) {
    const auto& parameter_types =
        MemoryManager::GetInstance()->getMemoryCache().parameter_types;

    if (std::find(parameter_types.begin(), parameter_types.end(), ParameterType::BSM)
        != parameter_types.end()) {
        poa_bsm = std::make_unique<ParamOptimizerAdapter>(
            std::vector<ParameterType>{ParameterType::BSM}
        );
    }
}


ParamOptimizerAdapter& StatParamOptimizerProxy::optimizer_for(
    const ParamId& pid
) {
    if (!pid.type.has_value()) {
        throw std::invalid_argument(
            "StatParamOptimizerProxy requires a typed ParamId for block '"
            + pid.block + "' and code '" + pid.code.to_string() + "'."
        );
    }

    if (pid.type.value() == ParameterType::BSM) {
        if (!poa_bsm) {
            throw std::logic_error(
                "A BSM parameter was requested, but the active HyperIso model "
                "does not provide a BSM parameter store."
            );
        }
        return *poa_bsm;
    }

    return poa_standard;
}


void StatParamOptimizerProxy::set_value(
    const ParamId& pid,
    scalar_t value
) {
    optimizer_for(pid).set_value(pid.block, pid.code, value);
}


void StatParamOptimizerProxy::set_param(
    const ParamId& pid,
    std::shared_ptr<Parameter> parameter
) {
    optimizer_for(pid).set_param(
        pid.block,
        pid.code,
        std::move(parameter)
    );
}


void StatParamOptimizerProxy::remove(const ParamId& pid) {
    optimizer_for(pid).remove(pid.block, pid.code);
}


void StatParamOptimizerProxy::commit(bool coalesce) {
    try {
        // Commit the model point first.  The following standard-parameter
        // commit then evaluates all dependent quantities at that BSM point.
        if (poa_bsm) {
            poa_bsm->commit(coalesce);
        }
        poa_standard.commit(coalesce);
    } catch (...) {
        if (poa_bsm) {
            poa_bsm->clear();
        }
        poa_standard.clear();
        throw;
    }
}


void StatParamOptimizerProxy::clear() {
    if (poa_bsm) {
        poa_bsm->clear();
    }
    poa_standard.clear();
}
