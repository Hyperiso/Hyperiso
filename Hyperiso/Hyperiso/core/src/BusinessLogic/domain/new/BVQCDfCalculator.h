#ifndef BVQDCFCALCULATOR_H
#define BVQDCFCALCULATOR_H

#include "Include.h"
#include "BVFFCalculator.h"
#include "BaseQCDfCalculator.h"

class BVQCDfCalculator : public BaseQCDfCalculator {
private:
    std::shared_ptr<BVFFCalculator> ff_calculator;
    complex_t G_perp_cache {};
    complex_t H_perp_cache {};
    complex_t H_2_cache {};
    double H_8_cache {0.0};
    bool static_integrals_ready {false};
    std::map<std::pair<bool, double>, complex_t> i_hsa1_cache;
    std::map<std::pair<bool, double>, complex_t> i_hsa2_cache;

    struct IndexedWilsonCache {
        complex_t C1{};
        complex_t C2{};
        complex_t C3{};
        complex_t C4{};
        complex_t C5{};
        complex_t C6{};
        complex_t C7{};
        complex_t C8{};
        complex_t C9{};
        complex_t C10{};
        complex_t CP7{};
        complex_t CP9{};
        complex_t CP10{};
        complex_t CQ1{};
        complex_t CQ2{};
        complex_t CPQ1{};
        complex_t CPQ2{};
    } wc;

    void build_indexed_wilsons(const std::map<WCoef, complex_t>& C);
    void precompute_static_integrals();
    double xi_perp_p_at(double q2);
    double xi_perp_m_at(double q2);

public:
    struct TTriplet {
        complex_t perp_p {};
        complex_t perp_m {};
        complex_t par_m {};
    };

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

    TTriplet evaluate_T_triplet(double q2, bool bar);

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


#endif // BVQDCFCALCULATOR_H
