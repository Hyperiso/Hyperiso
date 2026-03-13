#include "StatParamOptimizerProxy.h"

/**
 * @brief Constructs the proxy with the parameter types used by the statistics layer.
 *
 * The underlying @ref ParamOptimizerAdapter is initialized to operate on:
 * - SM parameters,
 * - FLAVOR parameters,
 * - DECAY parameters,
 * - WILSON parameters.
 */
StatParamOptimizerProxy::StatParamOptimizerProxy() : poa({ParameterType::SM, ParameterType::FLAVOR, ParameterType::DECAY, ParameterType::WILSON}) { //TODO : add the right parameterType
    
}

/**
 * @copydoc IStatParamOptimizerProxy::set_value
 */
void StatParamOptimizerProxy::set_value(const BlockName& block, const LhaID& id, scalar_t v) {
    poa.set_value(block, id, v);
}

/**
 * @copydoc IStatParamOptimizerProxy::set_param
 */
void StatParamOptimizerProxy::set_param(const BlockName& block, const LhaID& id, std::shared_ptr<Parameter> p) {
    poa.set_param(block,id, p);
}

/**
 * @copydoc IStatParamOptimizerProxy::remove
 */
void StatParamOptimizerProxy::remove(const BlockName& block, const LhaID& id) {
    poa.remove(block, id);
}

/**
 * @copydoc IStatParamOptimizerProxy::commit
 */
void StatParamOptimizerProxy::commit(bool coalesce) {
    poa.commit(coalesce);
}

/**
 * @copydoc IStatParamOptimizerProxy::clear
 */
void StatParamOptimizerProxy::clear() {
    poa.clear();
}