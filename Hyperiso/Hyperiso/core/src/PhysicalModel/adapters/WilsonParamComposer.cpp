#include "WilsonParamComposer.h"

#include <algorithm>
#include <vector>

void WilsonParamComposer::compose_block(
    const std::string &block_name,
    const std::unordered_map<ParameterType, std::vector<std::string>> &source_names,
    const DepUpdateFunc &update_func) {
    CompositeParamAdapter().add_block_dependency(block_name, source_names, ParameterType::WILSON, update_func);
    composed_blocks.emplace(block_name);
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
    composed_blocks.emplace(pid.block);
}

void WilsonParamComposer::remove_block(const std::string &block_name) {
    CompositeParamAdapter().remove_dependency(block_name, ParameterType::WILSON);
    composed_blocks.erase(block_name);
}

void WilsonParamComposer::update(const std::string &block_name) {
    CompositeParamAdapter().update_dependency(block_name, ParameterType::WILSON);
}

void WilsonParamComposer::remove_all_composed_blocks() {
    // The dependency graph is global in Core::Parameters, while this registry is
    // manager/composer-local.  Remove downstream blocks before their sources so
    // no lazy DependentBlock keeps a stale pointer to a source block erased from
    // the global repository.
    std::vector<std::string> blocks(composed_blocks.begin(), composed_blocks.end());

    auto removal_rank = [](const std::string& name) {
        // Final hadronic Wilson blocks depend on the intermediate running
        // blocks, which depend on matching/helper blocks. Remove in that order.
        if (name.find("_B_SCALE_") != std::string::npos &&
            name.find("__") == std::string::npos) {
            return 0;
        }
        if (name.find("__") != std::string::npos) {
            return 1;
        }
        if (name.find("_EW_SCALE") != std::string::npos) {
            return 2;
        }
        return 3;
    };

    std::sort(blocks.begin(), blocks.end(), [&](const std::string& a, const std::string& b) {
        const int ra = removal_rank(a);
        const int rb = removal_rank(b);
        if (ra != rb) return ra < rb;
        return a < b;
    });

    for (const auto& block_name : blocks) {
        CompositeParamAdapter().remove_dependency(block_name, ParameterType::WILSON);
    }
    composed_blocks.clear();
}
