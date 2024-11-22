#ifndef __B__KSTAR_GAMMA_H__
#define __B__KSTAR_GAMMA_H__

#include "Observable.h"

class Delta0_B__Kstar_gamma : public Observable {
public:
    Delta0_B__Kstar_gamma(Model model, QCDOrder order, double scale) : Observable(Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, 1.09e-10, 7.4e-11, model, order, scale) {};
    double eval() const override;
};

#endif // __B__KSTAR_GAMMA_H__