#ifndef RARE_DECAY_CALCULATOR_H
#define RARE_DECAY_CALCULATOR_H

#include <complex>
#include <memory>

#include "Observable.h"
#include "Wilson.h"

typedef std::complex<double> complex_t;

class ObsEvaluator {
private:
    static WilsonManager* computeWilsons(int model, int order, double scale);

    static complex_t Bs_mumu(WilsonManager* wm);
    static complex_t Bd_mumu(WilsonManager* wm);

public:
    static complex_t Evaluate(Observable* o);
    
};

#endif // RARE_DECAY_CALCULATOR_H
