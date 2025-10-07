#if !defined(HYPERISO_OBSERVABLE_H) 
#define HYPERISO_OBSERVABLE_H

#include "Include.h"
#include "Compound.h"
#include "DecayParent.h"

class Observable : public Compound {

protected:
    const ObservableId id;
    std::shared_ptr<DecayParent> decay_parent;

public:
    Observable(ObservableId id, std::shared_ptr<DecayParent> decay_parent) : id(id), decay_parent(decay_parent) {}

    ObservableId getId() const { return id; }
    scalar_t get_exp_val() const;
    scalar_t get_exp_uncertainty(UncertaintyType u_type=UncertaintyType::COMBINED) const;
    Estimate get_exp() const;
    scalar_t eval() const override;
    std::vector<ObservableValue> compute() const;
}; 


#endif // HYPERISO_OBSERVABLE_H
