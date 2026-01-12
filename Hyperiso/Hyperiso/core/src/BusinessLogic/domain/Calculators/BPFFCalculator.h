#ifndef __BPFFCALCULATOR_H__
#define __BPFFCALCULATOR_H__

#include "Include.h"
#include "ObsParameterProxy.h"
#include "IFFCalculator.h"

enum class BP_FF {F_PLUS, F_0, F_T, XI_P};
enum class BP_FF_Src {AS, GRvDV, GKvD_SR_LAT, GKvD_SR, FLAG24, HPQCD22};

class BPFFCalculator : public IFFCalculator<BP_FF> {
private:
    double m_B;
    double m_P;
    double t_p, t_m, t_0, z_0;
    double L_chi;
    std::map<BP_FF, std::array<double, 4>> alpha_ai {};
    std::map<BP_FF, double> m_R {};
    std::string src_block;

    static inline const std::map<LhaID, std::string> allowed_decays {
        {{511, 311}, "B_K"},
        {{521, 321}, "B_K"}
    };

    static inline const std::map<BP_FF_Src, size_t> sse_order {
        {BP_FF_Src::AS, 3},
        {BP_FF_Src::GRvDV, 2},
        {BP_FF_Src::GKvD_SR_LAT, 2},
        {BP_FF_Src::GKvD_SR, 2},
        {BP_FF_Src::FLAG24, 3},
        {BP_FF_Src::HPQCD22, 3}
    };

public:
    BPFFCalculator() = default;
    BPFFCalculator(int B_id, int P_id,
                        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm,
                        BP_FF_Src src=BP_FF_Src::GKvD_SR_LAT);

    complex_t z(double t, double t_p, double t_0);
    double E(double q2);
    double get(BP_FF a, double q2) override;

private:
    void load_FF_params(BP_FF_Src src);
    double pole(double q2, double m_R);
    double F_a(BP_FF a, double q2);

};

#endif // __BPFFCALCULATOR_H__
