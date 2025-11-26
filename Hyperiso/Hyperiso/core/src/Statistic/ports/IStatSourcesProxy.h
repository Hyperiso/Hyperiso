#ifndef ISTAT_SOURCES_PROXY_H
#define ISTAT_SOURCES_PROXY_H

#include <unordered_set>
#include "Include.h"

class IStatSourcesProxy {
public:
    virtual ~IStatSourcesProxy() = default;

    virtual std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const = 0;

};

#endif