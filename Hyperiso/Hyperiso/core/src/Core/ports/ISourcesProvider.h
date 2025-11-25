#include <unordered_set>
#include "Include.h"

class ISourcesProvider {
public:
    virtual ~ISourcesProvider() = default;

    virtual std::unordered_set<ParamId> get_all_leaf_sources(const std::unordered_set<ParamId>& param_ids) 
                                      const = 0;

};