#if !defined(HYPERISO_OBSERVABLE_H) 
#define HYPERISO_OBSERVABLE_H

#include "General.h"
#include "Compound.h"
#include "Math.h"
#include "WilsonManager.h"
#include "Parameters.h"
#include "epsilon_calculator.h"
#include "Wilson_THDMv2.h"
#include "Wilson_susyv2.h"
#include "WilsonManager.h"
#include "Wilsonv2.h"

class Observable : public Compound {

protected:

    const Observables id;
    const double exp_val;
    const double exp_std;
    Model model;
    QCDOrder order;
    double scale;

public:

    Observable(Observables id, double exp_val, double exp_std, Model model, QCDOrder order, double scale) 
        : id(id), exp_val(exp_val), exp_std(exp_std), model(model), order(order), scale(scale) {} 

    Observables getId() const;
    double get_exp_val() const;
    double get_exp_var() const;
    virtual double eval() const override = 0;
    CoefficientManager* computeWilsons(bool traditional_basis=false) const;
    CoefficientManager* computeWilsons(Model model, QCDOrder order, double scale, bool traditional_basis=false) const;

}; 


#endif // HYPERISO_OBSERVABLE_H
