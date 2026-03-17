#include "StatParameterProxy.h"

/**
 * @brief Constructs the statistics parameter proxy.
 *
 * Initializes:
 * - a generic provider (`pp`) for typed ParamId access,
 * - a typed provider (`pp_with_type`) bound to @p type for `(block, LhaID)` access.
 *
 * @param type Parameter namespace used by block-based queries.
 */
StatParameterProxy::StatParameterProxy(ParameterType type) { 
    if (!StatParameterProxy::ALLOWED.contains(type)) {
        LOG_ERROR("ValueError", "BusinessLogic cannot access parameter type", ParameterTypeMapper::str(type));
    }
    this->pp_with_type = ParameterProvider(type);
    this->pp = ParameterProvider();
}

/**
 * @copydoc IStatParameterProxy::get_param(const ParamId&) const
 */
std::shared_ptr<Parameter> StatParameterProxy::get_param(const ParamId& pid) const {
    return pp.get_parameter(pid);
}

/**
 * @copydoc IStatParameterProxy::get_param(const std::string&, const LhaID&) const
 */
std::shared_ptr<Parameter> StatParameterProxy::get_param(const std::string& block, const LhaID& id) const {
    return pp_with_type.get_parameter(ParamId(pp_with_type.get_type(), block, id));
}

/**
 * @copydoc IStatParameterProxy::operator()(const ParamId&, DataType) const
 */
scalar_t StatParameterProxy::operator()(const ParamId& pid, DataType d_type) const {
    if (!pid.type.has_value()) {
        LOG_WARN("LogicError", "Use of untyped ParamId in ParameterProvider.");
    }

    if (!StatParameterProxy::ALLOWED.contains(pid.type.value())) {
        LOG_ERROR("ValueError", "BusinessLogic cannot access parameter type", ParameterTypeMapper::str(pid.type.value()));
    }
    if (pid.type == ParameterType::WILSON) {
        return pp.exists(pid) ? pp(pid, d_type) : scalar_t();
    } 
    return pp(pid, d_type); 

}

/**
 * @copydoc IStatParameterProxy::operator()(const ObservableId&, DataType) const
 */
std::vector<double> StatParameterProxy::operator()(const ObservableId& id, DataType d_type) const {
    

    std::unordered_set<std::string> blocks = BlockProvider().get_all_blocks(ParameterType::OBSERVABLE);
    std::vector<double> out;
    for (auto block : blocks) {
        out.push_back(pp(ParamId(ParameterType::OBSERVABLE, block, ObservableMapper::flha(id)), d_type));
    }
    return out;
    // return pp(PKCaramId(ParameterType::OBSERVABLE, "FOBS", ObservableMapper::flha(id)), d_type);
}

/**
 * @copydoc IStatParameterProxy::operator()(const BinnedObservableId&, DataType) const
 */
std::vector<double> StatParameterProxy::operator()(const BinnedObservableId &id, DataType d_type) const {

    std::unordered_set<std::string> blocks = BlockProvider().get_all_blocks(ParameterType::OBSERVABLE);

    std::vector<double> out;
    for (auto block : blocks) {
        out.push_back(pp(ParamId(ParameterType::OBSERVABLE, block, id.flha()), d_type));
    }
    return out;
    // return pp(ParamId(ParameterType::OBSERVABLE, "FOBS", id.flha()), d_type);
}

/**
 * @copydoc IStatParameterProxy::operator()(const std::string&, const LhaID&, DataType) const
 */
scalar_t StatParameterProxy::operator()(const std::string& block, const LhaID& id, DataType d_type) const {
    if (pp_with_type.get_type() == ParameterType::WILSON) {
        return pp_with_type.exists(block, id) ? pp_with_type(block, id, d_type) : scalar_t();
    } 
    scalar_t value = pp_with_type(block, id, d_type);
    return value;
}

/**
 * @copydoc IStatParameterProxy::get_obs_param(const BinnedObservableId&) const
 */
std::vector<std::shared_ptr<Parameter>> StatParameterProxy::get_obs_param(const BinnedObservableId& id) const {

    std::unordered_set<std::string> blocks = BlockProvider().get_all_blocks(ParameterType::OBSERVABLE);

    std::vector<std::shared_ptr<Parameter>> out;
    for (auto block : blocks) {
        out.push_back(pp.get_parameter(ParamId(ParameterType::OBSERVABLE, block, id.flha())));
    }
    return out;
    // return pp.get_parameter(ParamId(ParameterType::OBSERVABLE, "FOBS", id.flha()));
}
