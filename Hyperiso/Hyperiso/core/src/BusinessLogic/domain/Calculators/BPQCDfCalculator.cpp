#include "BPQCDfCalculator.h"

BPQCDfCalculator::BPQCDfCalculator(int B_id, int V_id, double mu_b, const std::map<WCoef, complex_t> &C, std::shared_ptr<BPFFCalculator> ff_calculator, B_FF_Type ff_tp,
                        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p,
                        std::shared_ptr<IObsQCDProxy> iobs_qcdp) : 
    BaseQCDfCalculator(B_id, V_id, mu_b, C, ff_tp, p, iobs_qcdp)
{
    this->ff_calculator = ff_calculator;
}

complex_t BPQCDfCalculator::T_P(double q2, bool bar) {
    complex_t C_P;

    if (this->ff_tp == B_FF_Type::SOFT) {
        C_P = -(C_par_0(q2, 1, bar) + this->alpha_s_mu_b * iobs_qcdp->get_constants()->C_F / (4. * PI) * (C_par_f(q2, 1, bar) + C_par_nf(q2, bar)));
    } else {
        C_P = -this->alpha_s_mu_b * iobs_qcdp->get_constants()->C_F / (4. * PI) * C_par_nf(q2, bar);
    }

    // printf("C_par_nf = %.4e + %.4e i\n", real(C_par_nf(q2, bar)), imag(C_par_nf(q2, bar)));
    // printf("C_P = %.4e + %.4e i\n", real(C_P), imag(C_P));
    // printf("xi_P = %.4e\n", this->ff_calculator->get(BP_FF::XI_P, q2));
    // printf("pref_par = %.4e\n", this->pref_par);
    // printf("m_P = %.4e\n", this->m_X);
    // printf("I_P = %.4e + %.4e i\n", real(I_par_p(q2, bar)), imag(I_par_p(q2, bar)));
        
    return this->ff_calculator->get(BP_FF::XI_P, q2) * C_P - this->pref_par * I_par_p(q2, bar);
}

double BPQCDfCalculator::Delta_P_0(double q2) {
    return 1. + this->alpha_s_mu_f * iobs_qcdp->get_constants()->C_F / (2. * PI) * (
        1. - L(q2) + 3. * PI2 * (m_B - 2 * E(q2)) * this->f_B * this->f_X_par / (iobs_qcdp->get_constants()->Nc * this->lambda_B_p * this->ff_calculator->get(BP_FF::XI_P, q2) * std::pow(this->E(q2), 2)) * (1 + this->a_1_par + this->a_2_par)
    );
}

double BPQCDfCalculator::Delta_P_T(double q2) {
    return 1. + this->alpha_s_mu_b * iobs_qcdp->get_constants()->C_F / (2. * PI) * (
        L(q2) - this->L_b + 6. * PI2 * this->f_B * this->f_X_par / (iobs_qcdp->get_constants()->Nc * this->lambda_B_p * this->ff_calculator->get(BP_FF::XI_P, q2) * this->E(q2)) * (1 + this->a_1_par + this->a_2_par)
    );
}
