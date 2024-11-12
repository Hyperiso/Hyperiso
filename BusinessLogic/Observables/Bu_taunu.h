#include "Observable.h"

class BR_Bu_taunu : public Observable {
public:
    BR_Bu_taunu(int model, int order, double scale) : Observable(Observables::BR_BU_TAUNU, 1.09e-10, 7.4e-11, model, order, scale) {};
    double eval() const override;
};