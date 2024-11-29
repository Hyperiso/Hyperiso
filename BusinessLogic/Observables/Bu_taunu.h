#ifndef __BU_TAUNU_H__
#define __BU_TAUNU_H__

#include "Observable.h"

class BR_Bu_taunu : public Observable {
public:
    BR_Bu_taunu(Model model, QCDOrder order, double scale) : Observable(Observables::BR_BU_TAUNU, 1.09e-4, 2.4e-5, model, order, scale) {
        add_dependence({"MASS", 2});
        add_dependence({"MASS", 5});
        add_dependence({"RECKM", 02});
        add_dependence({"FMASS", 521});
        add_dependence({"FLIFE", 521});
        add_dependence({"FCONST", 52101});
    }
    double eval() const override;
};

#endif // __BU_TAUNU_H__
