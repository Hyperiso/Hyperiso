#ifndef RARE_DECAY_CALCULATOR_H
#define RARE_DECAY_CALCULATOR_H

#include <complex>
#include <memory>

#include "Observable.h"

typedef std::complex<double> complex_t;

class ObsEvaluator {
private:
    static double alpha_s;
    static complex_t Bs_mumu();
    static complex_t Bd_mumu();

public:
    static complex_t Evaluate(Observable* o);
};

#endif // RARE_DECAY_CALCULATOR_H
