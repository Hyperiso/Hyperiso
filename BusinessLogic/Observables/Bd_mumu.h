#ifndef __BD_MUMU_H__
#define __BD_MUMU_H__

#include "Observable.h"

class BR_Bd_mumu : public Observable {
public:
    BR_Bd_mumu(int model, int order, double scale) : Observable(Observables::BR_BD_MUMU, 1.09e-10, 7.4e-11, model, order, scale) {};
    double eval() const override;
};

#endif // __BD_MUMU_H__
