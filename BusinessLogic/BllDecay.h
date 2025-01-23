#ifndef __BLLDECAY_H__
#define __BLLDECAY_H__

#include "DecayParent.h"
#include "General.h"

/**
 * @brief Decay parent for the Bq > ll decays. Currently implements both CP-averaged and untagged Bs > mu+ mu- branching ratios and the CP-averaged Bd > mu+ mu- decays. 
 */
class BllDecay : public DecayParent {

protected:
    complex_t W1(double r, bool prime);
    complex_t W2Q(double r, bool prime);
    complex_t W210(double x, bool prime);
    double ckm(double V_tb_r, double V_tb_i, double V_ts_r, double V_ts_i);
    double BR_avg_Bq_mumu(complex_t w1, complex_t w2q, complex_t w210, double ckm, double b, double G_F, double inv_alpha, double f_Bs, double m_Bs, double life_Bs);
    double A_DG(double x, double r);
    double BR_untag_Bs_mumu(double br_avg, double ys, double A);

public:
    BllDecay(QCDOrder order, double matching_scale, double hadronic_scale) {
        winfo.matching_scale = matching_scale;
        winfo.hadronic_scale = hadronic_scale;
        winfo.model = MemoryManager::GetInstance()->getModel();
        winfo.order = order;
        winfo.basis = BWilsonBasis::STANDARD;
        winfo.wgroups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};

        build_op_tree();
    }

    void build_op_tree() override;

};

#endif // __BLLDECAY_H__