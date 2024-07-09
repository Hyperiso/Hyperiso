#if !defined(HYPERISO_OBSERVABLE_H) 
#define HYPERISO_OBSERVABLE_H

#include "Observables.h"
#include <complex>

typedef std::complex<double> complex_t;

class Observable {
private:
    Observables id;
    int order;
    int model;
    double scale;
    int wilson_basis;
    complex_t value;
    bool evaluated;

    void evaluate();  // Evaluates the observable at a given scale and sets its value, sets evaluated to false if evaluation fails

public:
    Observable(Observables id, double scale, int order, int model, int wilson_basis=1) : id(id), scale(scale), order(order), model(model), wilson_basis(wilson_basis) {
        this->evaluate();
    }

    complex_t getValue() const;  
    inline double getScale() const { return this->scale; }   
    inline Observables getId() const { return this->id; }   
    inline int getModel() const { return this->model; }   
    inline int getOrder() const { return this->order; }   
    inline int getWilsonBasis() const { return this->wilson_basis; }   
}; 


#endif // HYPERISO_OBSERVABLE_H
