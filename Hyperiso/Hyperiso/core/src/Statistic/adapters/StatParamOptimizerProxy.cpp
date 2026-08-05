#include "StatParamOptimizerProxy.h"

#include <stdexcept>
#include <utility>


StatParamOptimizerProxy::StatParamOptimizerProxy()
    : poa_bsm({ParameterType::BSM}),
      poa_standard({
          ParameterType::SM,
          ParameterType::FLAVOR,
          ParameterType::DECAY,
          ParameterType::WILSON
      }) {}


ParamOptimizerAdapter& StatParamOptimizerProxy::optimizer_for(
    const ParamId& pid
) {
    if (!pid.type.has_value()) {
        throw std::invalid_argument(
            "StatParamOptimizerProxy requires a typed ParamId for block '"
            + pid.block + "' and code '" + pid.code.to_string() + "'."
        );
    }

    return pid.type.value() == ParameterType::BSM
        ? poa_bsm
        : poa_standard;
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
        poa_bsm.commit(coalesce);
        poa_standard.commit(coalesce);
    } catch (...) {
        poa_bsm.clear();
        poa_standard.clear();
        throw;
    }
}


void StatParamOptimizerProxy::clear() {
    poa_bsm.clear();
    poa_standard.clear();
}
