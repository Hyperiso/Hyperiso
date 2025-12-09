#include "StatParamSourcesProxy.h"

std::unordered_set<ParamId> StatParamSourcesProxy::get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) const {
    return psp.get_all_leaf_sources(param_ids);
}