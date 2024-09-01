#ifndef RARE_DECAY_CALCULATOR_H
#define RARE_DECAY_CALCULATOR_H

#include <complex>
#include <memory>

#include "Observable.h"
#include "WilsonManager.h"

typedef std::complex<double> complex_t;

class ObsEvaluator {
private:
    void computeWilsons(int model, int order, double scale, bool traditional_basis);

    complex_t Bs_mumu(bool untag);
    complex_t Bd_mumu();
    complex_t Bu_taunu(int model, bool np_only);
    complex_t Delta_0_B_Kstargamma(double mu_b);

    CoefficientManager* manager;


public:
    complex_t Evaluate(Observable* o);
    
};

#endif // RARE_DECAY_CALCULATOR_H
