#if !defined(HYPERISO_OBSERVABLE_H) 
#define HYPERISO_OBSERVABLE_H

#include "Include.h"
#include "Compound.h"
#include "DecayParent.h"

class Observable : public Compound {

protected:
    const Observables id;
    std::shared_ptr<DecayParent> decay_parent;

public:
    Observable(Observables id, std::shared_ptr<DecayParent> decay_parent) : id(id), decay_parent(decay_parent) {}

    Observables getId() const { return id; }
    scalar_t get_exp_val() const;
    scalar_t get_exp_uncertainty(UncertaintyType u_type=UncertaintyType::COMBINED) const;
    Estimate get_exp() const;
    scalar_t eval() const override;
    size_t get_n_evals() const { return decay_parent->get_n_evals(id); };
}; 


#endif // HYPERISO_OBSERVABLE_H
