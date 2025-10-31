#ifndef __LBLFFCALCULATOR_H__
#define __LBLFFCALCULATOR_H__

#include "Include.h"
#include "ObsParameterProxy.h"
#include "IFFCalculator.h"

enum class LbL_FF {F_PERP, H_PERP, G_PERP, F_PLUS, H_PLUS, G_PLUS, H_TILDE_PERP, H_TILDE_PLUS};
enum class LbL_FF_Src {DM};

class LbLFFCalculator : public IFFCalculator<LbL_FF> {
private:
    double m_Lb, m_L;
    double t_p, t_0;
    std::map<LbL_FF, std::array<double, 3>> alpha_ai {};
    std::map<LbL_FF, double> m_R {};

    static inline const std::map<LbL_FF_Src, size_t> sse_order {
        {LbL_FF_Src::DM, 2},
    };

public:
    LbLFFCalculator(LbL_FF_Src src=LbL_FF_Src::DM);

    complex_t z(double t, double t_p, double t_0);
    double get(LbL_FF a, double q2) override;

private:
    void load_FF_params(LbL_FF_Src src);
    double pole(double q2, double m_R);
    double F_a(LbL_FF a, double q2);

};

#endif // __LBLFFCALCULATOR_H__
