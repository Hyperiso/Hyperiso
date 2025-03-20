#include "WilsonParamComposer.h"

void WilsonParamComposer::compose(
    const std::string &block_name,
    const std::unordered_map<ParameterType, std::vector<std::string>> &source_names,
    const DepUpdateFunc &update_func) 
{
    CompositeParamAdapter().add_dependency(block_name, source_names, ParameterType::WILSON, update_func);
}

void WilsonParamComposer::update(const std::string &block_name) {
    CompositeParamAdapter().update_dependency(block_name, ParameterType::WILSON);
}
