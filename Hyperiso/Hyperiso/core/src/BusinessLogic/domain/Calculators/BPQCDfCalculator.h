#ifndef __BPQCDFCALCULATOR_H__
#define __BPQCDFCALCULATOR_H__

#include "Include.h"
#include "BPFFCalculator.h"
#include "BaseQCDfCalculator.h"

class BPQCDfCalculator : public BaseQCDfCalculator {
private:
    std::shared_ptr<BPFFCalculator> ff_calculator;

public:
    BPQCDfCalculator() = default;
    BPQCDfCalculator(int B_id, int V_id, double mu_b, const std::map<WCoef, complex_t>& C, std::shared_ptr<BPFFCalculator> ff_calculator, B_FF_Type ff_tp,
                        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm,
                        std::shared_ptr<IObsQCDProxy> iobs_qcdp);

    complex_t T_P(double q2, bool bar);
    double Delta_P_0(double q2);
    double Delta_P_T(double q2);
};

#endif // __BPQCDFCALCULATOR_H__
