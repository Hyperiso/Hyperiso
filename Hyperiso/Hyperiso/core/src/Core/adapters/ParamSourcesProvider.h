#ifndef PARAM_SOURCES_PROVIDER_H
#define PARAM_SOURCES_PROVIDER_H

#include "ISourcesProvider.h"

class ParamSourcesProvider : public ISourcesProvider {
public:

    ParamSourcesProvider() = default;
    
    std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const override;

};

#endif