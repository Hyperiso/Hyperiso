#ifndef STAT_PARAM_SOURCES_PROXY_H
#define STAT_PARAM_SOURCES_PROXY_H

#include "IStatSourcesProxy.h"
#include "ParamSourcesProvider.h"

class StatParamSourcesProxy : public IStatSourcesProxy {
public:

    StatParamSourcesProxy() = default;
    
    std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const override;
private:
    ParamSourcesProvider psp;
};

#endif