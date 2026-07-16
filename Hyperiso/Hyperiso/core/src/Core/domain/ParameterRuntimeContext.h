#ifndef PARAMETER_RUNTIME_CONTEXT_H
#define PARAMETER_RUNTIME_CONTEXT_H

#include <map>
#include <memory>
#include <unordered_set>
#include <vector>

#include "BlockAccessor.h"
#include "config.hpp"

class Parameters;

/**
 * @brief Per-thread parameter runtime used to isolate Monte-Carlo workers.
 *
 * A context owns an independent input-cache clone and an independent cache of
 * Parameters repositories. When installed through ScopedParameterRuntimeContext,
 * Parameters::GetInstance routes all reads/writes to this local context instead
 * of the process-global singleton repositories.
 */
class ParameterRuntimeContext {
public:
    ParameterRuntimeContext();

    std::shared_ptr<Parameters> get_parameters(ParameterType id);
    std::shared_ptr<BlockAccessor> extract_blocks(const std::unordered_set<BlockName>& block_names) const;
    std::unordered_set<BlockName> available_input_blocks() const;
    const std::vector<ParameterType>& parameter_types() const;

    static ParameterRuntimeContext* current();

private:
    std::shared_ptr<BlockAccessor> local_input_cache_;
    std::map<ParameterType, std::shared_ptr<Parameters>> local_instances_;
    std::vector<ParameterType> parameter_types_;

    static thread_local ParameterRuntimeContext* current_;

    friend class ScopedParameterRuntimeContext;
};

/**
 * @brief RAII installer for a ParameterRuntimeContext on the current thread.
 */
class ScopedParameterRuntimeContext {
public:
    explicit ScopedParameterRuntimeContext(ParameterRuntimeContext& ctx);
    ~ScopedParameterRuntimeContext();

    ScopedParameterRuntimeContext(const ScopedParameterRuntimeContext&) = delete;
    ScopedParameterRuntimeContext& operator=(const ScopedParameterRuntimeContext&) = delete;

private:
    ParameterRuntimeContext* previous_ = nullptr;
};

#endif // PARAMETER_RUNTIME_CONTEXT_H
