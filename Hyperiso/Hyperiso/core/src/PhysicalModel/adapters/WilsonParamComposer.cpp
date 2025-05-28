#include "WilsonParamComposer.h"

void WilsonParamComposer::compose_block(
    const std::string &block_name,
    const std::unordered_map<ParameterType, std::vector<std::string>> &source_names,
    const DepUpdateFunc &update_func) {
    CompositeParamAdapter().add_block_dependency(block_name, source_names, ParameterType::WILSON, update_func);
    WilsonParamComposer::composed_blocks.emplace(block_name);
}

void WilsonParamComposer::compose_parameter(const ParamId& pid,
                                            const std::unordered_set<ParamId>& source_pids,
                                            const DepParamUpdateFunc& update_func)
{
    ParamId typed_pid = {ParameterType::WILSON, pid.block, pid.code};
    std::unordered_set<ParamId> typed_sources;
    for (auto& src_pid : source_pids) {
        if (!src_pid.type.has_value()){
            typed_sources.emplace(ParamId{ParameterType::WILSON, src_pid.block, src_pid.code});
        } else {
            typed_sources.emplace(ParamId{src_pid.type.value(), src_pid.block, src_pid.code});
        }
    }
    CompositeParamAdapter().add_param_dependency(typed_pid, typed_sources, update_func);
    WilsonParamComposer::composed_blocks.emplace(pid.block);
}

void WilsonParamComposer::remove_block(const std::string &block_name) {
    CompositeParamAdapter().remove_dependency(block_name, ParameterType::WILSON);
    WilsonParamComposer::composed_blocks.erase(block_name);
}

void WilsonParamComposer::update(const std::string &block_name) {
    CompositeParamAdapter().update_dependency(block_name, ParameterType::WILSON);
}

void WilsonParamComposer::remove_all_composed_blocks() {
    for (const auto& block_name : WilsonParamComposer::composed_blocks) {
        CompositeParamAdapter().remove_dependency(block_name, ParameterType::WILSON);
    }
    WilsonParamComposer::composed_blocks.clear();
}
