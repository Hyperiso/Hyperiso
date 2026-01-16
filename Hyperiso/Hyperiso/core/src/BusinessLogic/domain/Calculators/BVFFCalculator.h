#ifndef BVCOMMON_H
#define BVCOMMON_H

#include "Include.h"
#include "IObsParameterProxy.h"
#include "IFFCalculator.h"

enum class BV_FF {A0, A1, A12, V, T1, T2, T23, A2, T3, XI_PERP, XI_PAR, F_PERP, F_PAR, F_0};
enum class BV_FF_Src {BSZ_SR_LAT, BSZ_SR, GRvDV, GKvD_SR_LAT, GKvD_SR, HLMW};

class BVFFCalculator : public IFFCalculator<BV_FF> {
private:
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm;
    double m_B, m_B2, m_B4;
    double m_V, m_V2, m_V4;
    double t_p, t_m, t_0, z_0;
    std::map<BV_FF, std::array<double, 3>> alpha_ai;
    std::map<BV_FF, double> m_R;
    std::string src_block;

    static inline const std::map<LhaID, std::string> allowed_decays {
        {{511, 313}, "B_Ks"},
        {{521, 323}, "B_Ks"},
        {{531, 333}, "B_phi"}
    };

public:
    BVFFCalculator() = default;
    BVFFCalculator(int B_id, int V_id, std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm, BV_FF_Src src=BV_FF_Src::BSZ_SR_LAT);

    complex_t z(double t, double t_p, double t_0);
    double E(double q2);
    double get(BV_FF a, double q2) override;

private:
    void load_FF_params(BV_FF_Src src);
    double pole(double q2, double m_R);
    double F_a(BV_FF a, double q2);
    double A_2(double q2);
    double T_3(double q2);
    double xi_perp(double q2);
    double xi_par(double q2);
    double f_perp(double q2);
    double f_par(double q2);
    double f_0(double q2);
};

#endif // BVCOMMON_H
