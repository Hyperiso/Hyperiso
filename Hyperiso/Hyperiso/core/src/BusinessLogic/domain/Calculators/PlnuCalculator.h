#ifndef __PLNUCALCULATOR_H__
#define __PLNUCALCULATOR_H__

#include "Include.h"
#include "ObsParameterProxy.h"

class PlnuCalculator {
private:
    double m_P, tau_P, f_P;
    double m_l, m_qu, m_qd;
    double V_sq;
    complex_t C_A, C_P;

    static inline const std::map<int, std::array<int, 4>> allowed_P {
    //    P    u  d  ckm_u ckm_d
        {521, {2, 5, 0,    2}},
        {411, {4, 1, 1,    0}},
        {431, {4, 3, 1,    1}},
        {321, {2, 3, 0,    1}},
        {211, {2, 1, 0,    0}}
    };

public:
    PlnuCalculator() = default;
    PlnuCalculator(int P_id, int l_id, complex_t C_A, complex_t C_P,
        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm);

    double BR_0_SM();
    double R_SM_BSM();
};

#endif // __PLNUCALCULATOR_H__
