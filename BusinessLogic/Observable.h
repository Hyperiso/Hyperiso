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
    complex_t value;
    bool evaluated;

    void evaluate();  // Evaluates the observable at a given scale and sets its value, sets evaluated to false if evaluation fails

public:
    Observable(Observables id, double scale, int order, int model) : id(id), scale(scale), order(order), model(model) {
        this->evaluate();
    }

    complex_t getValue() const;  
    inline double getScale() const { return this->scale; }   
    inline Observables getId() const { return this->id; }   
    inline int getModel() const { return this->model; }   
    inline int getOrder() const { return this->order; }   
}; 


#endif // HYPERISO_OBSERVABLE_H
