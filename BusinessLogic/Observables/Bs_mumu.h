#ifndef __BS_MUMU_H__
#define __BS_MUMU_H__

#include "Observable.h"

class BR_Bs_mumu : public Observable {
public:
    BR_Bs_mumu(Model model, QCDOrder order, double scale) : Observable(Observables::BR_BS_MUMU, 2.4e-9, 4e-10, model, order, scale) {
        add_dependence({"MASS", 3});
        add_dependence({"MASS", 5});
        add_dependence({"MASS", 6});
        add_dependence({"RECKM", 22});
        add_dependence({"RECKM", 21});
        add_dependence({"FMASS", 531});
        add_dependence({"FLIFE", 531});
        add_dependence({"FCONST", 53101});
    };
    double eval() const override;
};

class BR_Bs_mumu_untag : public Observable {
public:
    BR_Bs_mumu_untag(Model model, QCDOrder order, double scale) : Observable(Observables::BR_BS_MUMU_UNTAG, 2.65e-9, 4.3e-10, model, order, scale) {};
    double eval() const override;
};

#endif // __BS_MUMU_H__