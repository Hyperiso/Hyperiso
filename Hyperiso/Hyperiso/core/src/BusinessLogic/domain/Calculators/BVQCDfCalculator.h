#ifndef __BVQDCFCALCULATOR_H__
#define __BVQDCFCALCULATOR_H__

#include "Include.h"
#include "BVFFCalculator.h"
#include "BaseQCDfCalculator.h"

class BVQCDfCalculator : public BaseQCDfCalculator {
private:
    std::shared_ptr<BVFFCalculator> ff_calculator;

public:
    BVQCDfCalculator() = default;
    BVQCDfCalculator(int B_id, int V_id, double mu_b, const std::map<WCoef, complex_t>& C, std::shared_ptr<BVFFCalculator> ff_calculator, B_FF_Type ff_tp,
                        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm,
                        std::shared_ptr<IObsQCDProxy> iobs_qcdp);

    double F_perp(double s_hat);
    double X_perp(double s_hat);
    complex_t H_perp();
    complex_t G_perp();
    complex_t H_2();
    double H_8();

    complex_t T_perp_p(double q2, bool bar);
    complex_t T_perp_m(double q2, bool bar);
    complex_t T_par_p(double q2, bool bar);
    complex_t T_par_m(double q2, bool bar);
    complex_t Delta_par(double q2);

    complex_t I_HSA_1(double q2, bool bar);
    complex_t I_HSA_2(double q2, bool bar);
    complex_t delta_T_perp_WA(double q2, bool bar);
    complex_t delta_T_perp_HSA(double q2, bool bar);
};


#endif // __BVQDCFCALCULATOR_H__
