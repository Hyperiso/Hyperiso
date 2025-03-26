#include "WilsonParamComposer.h"

void WilsonParamComposer::compose_block(
    const std::string &block_name,
    const std::unordered_map<ParameterType, std::vector<std::string>> &source_names,
    const DepUpdateFunc &update_func) 
{
    CompositeParamAdapter().add_block_dependency(block_name, source_names, ParameterType::WILSON, update_func);
}

void WilsonParamComposer::compose_parameter(const ParamId& pid,
                                            const std::unordered_set<ParamId>& source_pids,
                                            const DepParamUpdateFunc& update_func)
{
    CompositeParamAdapter().add_param_dependency(pid, source_pids, update_func);
}

void WilsonParamComposer::update(const std::string &block_name) {
    CompositeParamAdapter().update_dependency(block_name, ParameterType::WILSON);
}
