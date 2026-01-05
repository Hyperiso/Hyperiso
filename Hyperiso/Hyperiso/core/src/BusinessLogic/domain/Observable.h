#if !defined(HYPERISO_OBSERVABLE_H) 
#define HYPERISO_OBSERVABLE_H

#include "Include.h"
#include "DecayParent.h"

class Observable {

protected:
    const ObservableId id;
    std::shared_ptr<DecayParent> decay_parent;
    std::unordered_set<ParamId> dependences;

public:
    Observable(ObservableId id, std::shared_ptr<DecayParent> decay_parent) : id(id), decay_parent(decay_parent) {}

    ObservableId getId() const { return id; }
    scalar_t get_exp_val() const;
    scalar_t get_exp_uncertainty(UncertaintyType u_type=UncertaintyType::COMBINED) const;
    std::vector<ObservableValue> compute() const;

    void add_dependence(const ParamId& param_name);
    void add_dependences(const std::unordered_set<ParamId>& param_names);
    const std::unordered_set<ParamId>& get_dependences() const;
}; 


#endif // HYPERISO_OBSERVABLE_H
