#if !defined(HYPERISO_OBSERVABLE_H) 
#define HYPERISO_OBSERVABLE_H

#include "General.h"
#include "Compound.h"
#include "DecayParent.h"

class Observable : public Compound {

protected:
    const Observables id;
    double exp_val;
    std::shared_ptr<DecayParent> decay_parent;

public:
    Observable(Observables id, std::shared_ptr<DecayParent> decay_parent, double exp_val = 0.0) : id(id), decay_parent(decay_parent), exp_val(exp_val) {}

    Observables getId() const { return id; }
    double get_exp_val() const { return exp_val; }
    void set_exp_val(double value) { exp_val = value; }
    scalar_t eval() const override;
    size_t get_n_evals() const { return decay_parent->get_n_evals(id); };
}; 


#endif // HYPERISO_OBSERVABLE_H
