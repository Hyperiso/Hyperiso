#ifndef BNUNU_WILSON_H
#define BNUNU_WILSON_H

#include "Wilson.h"

/**
 * @brief WET Wilson coefficient for b -> s nu_i anti-nu_j.
 *
 * Convention follows arXiv:2301.06990:
 *   L_eff = 4 G_F/sqrt(2) lambda_t sum C_a O_a,
 *   O_{L,R}^{ij} = e^2/(4 pi)^2
 *       (sbar gamma_mu P_{L,R} b)
 *       (nubar_i gamma^mu (1-gamma5) nu_j).
 *
 * The phenomenological SM value C_L^SM=-6.32 already includes the quoted
 * NLO-QCD and two-loop electroweak corrections.  It is therefore stored as
 * the LO/base coefficient in HyperIso's additive perturbative convention;
 * the NLO and NNLO increments are zero.
 */
class BNuNuWilson : public WilsonCoefficient {
public:
    BNuNuWilson(WCoef coefficient, ContributionType contribution);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<BNuNuWilson>(*this);
    }

    static bool is_left(WCoef coefficient);
    static bool is_sm_diagonal_left(WCoef coefficient);

private:
    WCoef coefficient_;
};

#endif
