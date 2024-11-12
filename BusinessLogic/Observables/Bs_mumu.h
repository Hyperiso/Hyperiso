#include "Observable.h"

class BR_Bs_mumu : public Observable {
public:
    BR_Bs_mumu(int model, int order, double scale) : Observable(Observables::BR_BS_MUMU, 2.4e-9, 4e-10, model, order, scale) {};
    double eval() const override;
};

class BR_Bs_mumu_untag : public Observable {
public:
    BR_Bs_mumu_untag(int model, int order, double scale) : Observable(Observables::BR_BS_MUMU_UNTAG, 2.65e-9, 4.3e-10, model, order, scale) {};
    double eval() const override;
};