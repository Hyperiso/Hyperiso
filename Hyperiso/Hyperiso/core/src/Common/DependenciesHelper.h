#pragma once
#include "Include.h"

class DependenciesHelper {
public:
    static std::unordered_set<ParamId> get_allowed_parameters(Observables id);
    static bool is_param_allowed(Observables id, ParamId pid);

    static std::unordered_set<ParamId> get_allowed_parameters(ObservableId id);
    static bool is_param_allowed(ObservableId id, ParamId pid);
    

private:
    static const std::map<DecayId, std::unordered_set<ParamId>> dep_lists;
};
