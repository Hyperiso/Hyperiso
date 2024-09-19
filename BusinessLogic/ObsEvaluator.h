#ifndef RARE_DECAY_CALCULATOR_H
#define RARE_DECAY_CALCULATOR_H

#include <complex>
#include <memory>

#include "Observable.h"
#include "WilsonManager.h"

typedef std::complex<double> complex_t;

class ObsEvaluator {
private:
    static CoefficientManager* computeWilsons(int model, int order, double scale, bool traditional_basis);

    static complex_t Bs_mumu(CoefficientManager* manager, bool untag);
    static complex_t Bd_mumu(CoefficientManager* manager);
    static complex_t Bu_taunu(int model, bool np_only);
    static complex_t Delta_0_B_Kstargamma(CoefficientManager* manager, double mu_b);


public:
    static complex_t Evaluate(Observable* o);
    
};

#endif // RARE_DECAY_CALCULATOR_H
