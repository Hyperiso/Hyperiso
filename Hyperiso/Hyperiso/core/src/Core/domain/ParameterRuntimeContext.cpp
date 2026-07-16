#include "ParameterRuntimeContext.h"

#include "MemoryManager.h"
#include "Parameters.h"

#include <mutex>

thread_local ParameterRuntimeContext* ParameterRuntimeContext::current_ = nullptr;

namespace {
std::mutex parameter_runtime_snapshot_mutex;
}

ParameterRuntimeContext::ParameterRuntimeContext() {
    std::lock_guard<std::mutex> lock(parameter_runtime_snapshot_mutex);
    auto* mm = MemoryManager::GetInstance();
    const auto& cache = mm->getMemoryCache();
    parameter_types_ = cache.parameter_types;
    local_input_cache_ = mm->clone_input_cache_deep();

    // Preserve runtime edits made on global plain input blocks before entering
    // the worker context, but do not copy dependent graph nodes. Dependent
    // blocks are rebuilt locally by the Parameters strategies.
    for (auto type : parameter_types_) {
        auto global_params = ParametersFactory::GetParameters(type);
        if (!global_params) {
            continue;
        }

        for (const auto& block_name : global_params->get_blocks_list()) {
            if (global_params->is_dependent_block(block_name)) {
                continue;
            }

            auto global_blocks = global_params->get_block_accessor();
            if (!global_blocks || !global_blocks->contains(block_name)) {
                continue;
            }

            auto cloned_block = global_blocks->at(block_name)->deep_clone_plain();
            local_input_cache_->emplace(block_name, cloned_block);
        }
    }
}

std::shared_ptr<Parameters> ParameterRuntimeContext::get_parameters(ParameterType id) {
    auto it = local_instances_.find(id);
    if (it != local_instances_.end()) {
        return it->second;
    }

    // The repository must be registered in the local map before
    // postInitialization runs. Post-initialization attaches dependent blocks via
    // DependentBlockManager, and that path can recursively call
    // Parameters::GetInstance(id) for the repository currently being built.
    // Registering only after CreateUncached returns causes unbounded recursion
    // and typically shows up as a segfault/stack overflow before MC accepts its
    // first draw.
    try {
        return ParametersFactory::CreateUncachedRegistered(
            id,
            [this, id](const std::shared_ptr<Parameters>& params) {
                local_instances_[id] = params;
            }
        );
    } catch (...) {
        local_instances_.erase(id);
        throw;
    }
}

std::shared_ptr<BlockAccessor> ParameterRuntimeContext::extract_blocks(
    const std::unordered_set<BlockName>& block_names
) const {
    return (*local_input_cache_)[block_names];
}

std::unordered_set<BlockName> ParameterRuntimeContext::available_input_blocks() const {
    return local_input_cache_->get_block_names();
}

const std::vector<ParameterType>& ParameterRuntimeContext::parameter_types() const {
    return parameter_types_;
}

ParameterRuntimeContext* ParameterRuntimeContext::current() {
    return current_;
}

ScopedParameterRuntimeContext::ScopedParameterRuntimeContext(ParameterRuntimeContext& ctx)
    : previous_(ParameterRuntimeContext::current_) {
    ParameterRuntimeContext::current_ = &ctx;
}

ScopedParameterRuntimeContext::~ScopedParameterRuntimeContext() {
    ParameterRuntimeContext::current_ = previous_;
}
