#include "BaseQCDfCalculator.h"

double BaseQCDfCalculator::E(double q2) {
    return (m_B * m_B + m_X * m_X - q2) / (2 * m_B);
}

BaseQCDfCalculator::BaseQCDfCalculator(int B_id, int X_id, double mu_b, const std::map<WCoef, complex_t> &C, B_FF_Type ff_tp,
                        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p,
                        std::shared_ptr<IObsQCDProxy> iobs_qcdp) :
    mu_b(mu_b), C(C), ff_tp(ff_tp), iobs_qcdp(iobs_qcdp)
{
    if (!this->allowed_decays.contains({B_id, X_id})) {
        LOG_ERROR("ValueError", "Wrong meson PDG code in BaseQCDfCalculator constructor:", B_id, ",", X_id);
    }

    this->src_block = allowed_decays.at({B_id, X_id});
    this->delta_qu = double (B_id == 521);
    this->fill_wilson_bar_cache();

    double beta_0 = iobs_qcdp->get_constants()->beta[4][0];
    auto run = [this, beta_0] (double value_1gev, double eta, double gamma) { return value_1gev * pow(eta, gamma / beta_0); };

    this->Lambda_h = (*p)(ParamId{ParameterType::DECAY, this->src_block, 14}, DataType::VALUE);
    double mu_f = sqrt(this->mu_b * this->Lambda_h);
    this->alpha_s_mu_b = (*iobs_qcdp)(AlphasConfig(this->mu_b, MassType::POLE, MassType::POLE));
    this->alpha_s_mu_f = (*iobs_qcdp)(AlphasConfig(mu_f, MassType::POLE, MassType::POLE));
    this->loop_f_mu_f = this->alpha_s_mu_f * iobs_qcdp->get_constants()->C_F / (4 * PI);
    this->loop_f_mu_b = this->alpha_s_mu_b * iobs_qcdp->get_constants()->C_F / (4 * PI);
    this->m_c_pole = (*p)(ParamId{ParameterType::SM, "QCD", {4, 1}}, DataType::VALUE);
    this->m_b_pole = (*p)(ParamId{ParameterType::SM, "QCD", {5, 5}}, DataType::VALUE);
    double eta_f = this->alpha_s_mu_f / (*iobs_qcdp)(AlphasConfig(1.0, MassType::POLE, MassType::POLE));
    double m_b_pole_2loop = (*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE);
    this->m_b_PS = m_b_pole_2loop - 4 * (*iobs_qcdp)(AlphasConfig(m_b_pole_2loop, MassType::POLE, MassType::POLE)) * mu_f / (3 * PI);
    this->m_B = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", B_id}, DataType::VALUE);
    this->m_Bd = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 511}, DataType::VALUE);
    this->m_X = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", X_id}, DataType::VALUE);
    this->f_B = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {B_id, 1}}, DataType::VALUE);
    this->f_X_par = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {X_id, 1}}, DataType::VALUE);
    this->lambda_B_p = (*p)(ParamId{ParameterType::DECAY, this->src_block, 13}, DataType::VALUE) / (1. - this->alpha_s_mu_f * log(pow(mu_f, 2)) * 1.8 / (3. * PI));
    this->lambda_hat_u = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE) 
                            / (std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE));
    // this->lambda_hat_u = 0.0; // ASK : Why neglected in SI for B > K l l ?
    this->a_1_par = run((*p)(ParamId{ParameterType::DECAY, this->src_block, {8, 1}}, DataType::VALUE), eta_f, gamma_par(1));
    this->a_2_par = run((*p)(ParamId{ParameterType::DECAY, this->src_block, {8, 2}}, DataType::VALUE), eta_f, gamma_par(2));

    // printf("a1par = %.4e\n", a_1_par);
    // printf("a2par = %.4e\n", a_2_par);

    this->e_q = B_id == 521 ? e_u : e_d;
    this->z_c = std::pow(this->m_c_pole / this->m_b_PS, 2);
    this->L_b = std::log(this->mu_b / this->m_b_PS);
    this->Delta_M = -6. * this->L_b - 4. * (1. - mu_f / this->m_b_PS);
    int Nc = iobs_qcdp->get_constants()->Nc;
    this->pref_par = PI2 * this->f_B * this->f_X_par / (Nc * this->m_B);

    if (X_id == 333) {
        this->n_T_par_m_0 = -4. * this->m_B / this->m_b_PS * (this->C_bar[WCoef::C3] + 3. * this->C_bar[WCoef::C4] + 12. * (this->C[WCoef::C3] + 10. * this->C[WCoef::C5]) - this->lambda_hat_u * (4. / 3. * this->C[WCoef::C1] + this->C[WCoef::C2]));
        this->n_T_par_m_0_bar = -4. * this->m_B / this->m_b_PS * (this->C_bar[WCoef::C3] + 3. * this->C_bar[WCoef::C4] + 12. * (this->C[WCoef::C3] + 10. * this->C[WCoef::C5]) - std::conj(this->lambda_hat_u) * (4. / 3. * this->C[WCoef::C1] + this->C[WCoef::C2]));
    } else {
        this->n_T_par_m_0 = 4. * this->m_B / this->m_b_PS * (3. * this->lambda_hat_u * this->delta_qu * this->C[WCoef::C2] - this->C_bar[WCoef::C3] - 3. * this->C_bar[WCoef::C4]);
        this->n_T_par_m_0_bar = 4. * this->m_B / this->m_b_PS * (3. * std::conj(this->lambda_hat_u) * this->delta_qu * this->C[WCoef::C2] - this->C_bar[WCoef::C3] - 3. * this->C_bar[WCoef::C4]);
    }

    // printf("n_T_par_m_0 = %.4e + %.4e i\n", real(e_q * n_T_par_m_0), imag(e_q * n_T_par_m_0));
    // printf("z_c = %.4e\n", z_c);

    bool isV = (X_id == 313 || X_id == 323 || X_id == 333);
    if (isV) {
        // printf("gamma1 = %.4e\n", gamma_perp(1));
        // printf("gamma2 = %.4e\n", gamma_perp(2));
        // printf("a10 = %.4e\n", std::real((*p)(ParamId{ParameterType::DECAY, this->src_block, {7, 1}})));
        // printf("a20 = %.4e\n", std::real((*p)(ParamId{ParameterType::DECAY, this->src_block, {7, 2}})));
        this->f_X_perp = run((*p)(ParamId{ParameterType::FLAVOR, "FCONST", {X_id, 2}}, DataType::VALUE), eta_f, iobs_qcdp->get_constants()->C_F);
        this->a_1_perp = run((*p)(ParamId{ParameterType::DECAY, this->src_block, {7, 1}}, DataType::VALUE), eta_f, gamma_perp(1));
        this->a_2_perp = run((*p)(ParamId{ParameterType::DECAY, this->src_block, {7, 2}}, DataType::VALUE), eta_f, gamma_perp(2));

        
        // printf("a_1_perp = %.4e\n", a_1_perp);
        // printf("a_2_perp = %.4e\n", a_2_perp);
        this->zeta_3_A = (*p)(ParamId{ParameterType::DECAY, this->src_block, 9}, DataType::VALUE);
        this->zeta_3_V = (*p)(ParamId{ParameterType::DECAY, this->src_block, 10}, DataType::VALUE);
        this->omega_10_A = (*p)(ParamId{ParameterType::DECAY, this->src_block, 11}, DataType::VALUE);
        this->delta_t_p = (*p)(ParamId{ParameterType::DECAY, this->src_block, {12, 1}}, DataType::VALUE);
        this->delta_t_m = (*p)(ParamId{ParameterType::DECAY, this->src_block, {12, 2}}, DataType::VALUE);
        this->pref_perp = PI2 * this->f_B * this->f_X_perp / (Nc * this->m_B);
    }

    // printf("m_b_PS = %.4e\n", m_b_PS);
    // printf("m_B = %.4e\n", m_B);

    // printf("lambda_Bp = %.4e\n", lambda_B_p);
    // printf("lambda_Bm = %.4e + %.4e i\n", real(1. / inv_lambda_B_m(1.0)), imag(1. / inv_lambda_B_m(1.0)));
    // printf("lambda_hat_u = %.4e + %.4e i\n", real(lambda_hat_u), imag(lambda_hat_u));

    // printf("C1bar = %.4e + %.4e i\n", real(C_bar[WCoef::C1]), imag(C_bar[WCoef::C1]));
    // printf("C2bar = %.4e + %.4e i\n", real(C_bar[WCoef::C2]), imag(C_bar[WCoef::C2]));
    // printf("C3bar = %.4e + %.4e i\n", real(C_bar[WCoef::C3]), imag(C_bar[WCoef::C3]));
    // printf("C4bar = %.4e + %.4e i\n", real(C_bar[WCoef::C4]), imag(C_bar[WCoef::C4]));
    // printf("C5bar = %.4e + %.4e i\n", real(C_bar[WCoef::C5]), imag(C_bar[WCoef::C5]));
    // printf("C6bar = %.4e + %.4e i\n", real(C_bar[WCoef::C6]), imag(C_bar[WCoef::C6]));
    
    // printf("phi_Kstar_perp = %.4e\n", phi_X(0.5, a_1_perp, a_2_perp));
    // printf("t_perp_mc(s = 1.0, u = 0.5) = %.4e + %.4e i\n", real(t_perp(0.5, this->m_c_pole, 1.0, E(1.0))), imag(t_perp(0.5, this->m_c_pole, 1.0, E(1.0))));
    // printf("t_perp_mb(s = 1.0, u = 0.5) = %.4e + %.4e i\n", real(t_perp(0.5, this->m_b_PS, 1.0, E(1.0))), imag(t_perp(0.5, this->m_b_PS, 1.0, E(1.0))));
    // printf("t_perp_0(s = 1.0, u = 0.5) = %.4e + %.4e i\n", real(t_perp(0.5, 0., 1.0, E(1.0))), imag(t_perp(0.5, 0., 1.0, E(1.0))));
    // printf("T_perp_p_nf(s = 1.0, u = 0.5) = %.4e + %.4e i\n", real(T_perp_p_nf(0.5, 1.0, false)), imag(T_perp_p_nf(0.5, 1.0, false)));

    // double v = this->m_B * this->m_B * (1 - 0.5) + 1.0 * 0.5;
    // printf("h_mc = %.4e + %.4e i\n", real(BV::h(v, this->m_c_pole, this->mu_b)), imag(BV::h(v, this->m_c_pole, this->mu_b)));
    // printf("h_mb(mb = %.4e) = %.4e + %.4e i\n", this->m_b_pole, real(BV::h(v, this->m_b_pole, this->mu_b)), imag(BV::h(v, this->m_b_pole, this->mu_b)));
    // printf("h_0 = %.4e + %.4e i\n", real(BV::h(v, 0, this->mu_b)), imag(BV::h(v, 0, this->mu_b)));
    // printf("T_par_m_nf = %.4e + %.4e i\n", real(e_q * T_par_m_nf(0.5, 1.0, false)), imag(e_q * T_par_m_nf(0.5, 1.0, false)));
    // printf("T_par_p_nf = %.4e + %.4e i\n", real(T_par_p_nf(0.5, 1.0, false)), imag(T_par_p_nf(0.5, 1.0, false)));

    // printf("T_par_m_0 = %.4e + %.4e i\n", real(T_par_m_0(false)), imag(T_par_m_0(false)));

}

double BaseQCDfCalculator::phi_X(double u, double a1, double a2) {
    double x = 2 * u - 1;
	double C1 = 3 * x;
	double C2 = -1.5 + 7.5 * x * x;

	return 6 * u * (1 - u) * (1 + a1 * C1 + a2 * C2); 
}

complex_t BaseQCDfCalculator::Y(double q2) {
    return BV::h(q2, this->m_c_pole, this->mu_b) * (4./3. * this->C[WCoef::C1] + this->C[WCoef::C2] + 6. * this->C[WCoef::C3] + 60. * this->C[WCoef::C5])
            - 0.5 * BV::h(q2, this->m_b_pole, this->mu_b) * (7. * this->C[WCoef::C3] + 4./3. * this->C[WCoef::C4] + 76. * this->C[WCoef::C5] + 64./3. * this->C[WCoef::C6])
            - 0.5 * BV::h(q2, 0., this->mu_b) * (this->C[WCoef::C3] + 4./3. * this->C[WCoef::C4] + 16. * this->C[WCoef::C5] + 64./3. * this->C[WCoef::C6])
            + 4./3. * this->C[WCoef::C3] + 64./9. * this->C[WCoef::C5] + 64./27. * this->C[WCoef::C6];
}

complex_t BaseQCDfCalculator::Y_u(double q2) {
    return (BV::h(q2, this->m_c_pole, this->mu_b) - BV::h(q2, 0., this->mu_b)) * (4. / 3. * this->C[WCoef::C1] + this->C[WCoef::C2]);
}

complex_t BaseQCDfCalculator::t_perp(double u, double m_q, double q2, double E_Kstar) {
    double mB2 = this->m_Bd * this->m_Bd;
    // ASK : Why always m_Bd and not m_B for B+ decay ? 
    if(fpeq(q2, 0.)) {
		if (fpeq(m_q, 0.)) return 4./(1.-u);
		double epsilon=1.e-10;
		complex_t xp=0.5+std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*mB2));
		complex_t xm=0.5-std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*mB2));
		return 4./(1.-u)*(1.+2.*m_q*m_q/(1.-u)/mB2*(BV::L_1(xp)+BV::L_1(xm)));
	} else {
        double s_hat = q2 / mB2;
        double mq_hat = m_q / m_Bd;
        return 2.*m_Bd/(1.-u)/E_Kstar*BV::I_1(u,s_hat,mq_hat)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(BV::B_0((1.-u)*mB2+u*q2,m_q)-BV::B_0(q2,m_q));
    }
}

complex_t BaseQCDfCalculator::t_par(double u, double m_q, double q2, double E_Kstar) {
    double mB2 = this->m_Bd * this->m_Bd;
    double s_hat = q2 / mB2;
    double mq_hat = m_q / this->m_Bd;
    
    // printf("In t_par(u = %.4e, m_q = %.4e, q2 = %.4e, E_K* = %.4e)\n", u, m_q, q2, E_Kstar);
    // printf("I_1 = %.4e + %.4e i\n", real(BV::I_1(u,s_hat,mq_hat)), imag(BV::I_1(u,s_hat,mq_hat)));
    // printf("B_0(u) = %.4e + %.4e i\n", real(BV::B_0((1.-u)*mB2+u*q2,m_q)), imag(BV::B_0((1.-u)*mB2+u*q2,m_q)));
    // printf("B_0(1) = %.4e + %.4e i\n", real(BV::B_0(q2,m_q)), imag(BV::B_0(q2,m_q)));

    // return 2.*this->m_B/(1.-u)/E_Kstar*BV::I_1(u,s_hat,mq_hat)+((1.-u)*mB2+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(BV::B_0((1.-u)*mB2+u*q2,m_q)-BV::B_0(q2,m_q));
    return 2.*this->m_Bd/(1.-u)/E_Kstar*BV::I_1(u,s_hat,mq_hat)+((1.-u)*mB2+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(BV::B_0((1.-u)*mB2+u*q2,m_q)-BV::B_0(q2,m_q));

}

complex_t BaseQCDfCalculator::T_par_m_0(bool bar) {
    return bar ? this->n_T_par_m_0_bar : this->n_T_par_m_0;
}

void BaseQCDfCalculator::fill_wilson_bar_cache() {
    auto b_ids = WCoefMapper::B_group();
    for (size_t i = 0; i < 6; i++) {
        this->C_bar[b_ids[i]] = 0;
        for (size_t j = 0; j < 6; j++) {
            this->C_bar[b_ids[i]] += P_bar[i][j] * this->C[b_ids[j]];
        }
    }
}

double BaseQCDfCalculator::gv_dga_4(double u) {
    double a1 = -60. * this->zeta_3_A * (this->omega_10_A + 4.) + 1680. * this->zeta_3_V;
    double a2 = 30. * this->zeta_3_A * (15. * this->omega_10_A + 32.) - 12600. * this->zeta_3_V + 36. * this->a_1_par - 72. * this->a_2_par - 12.;
    double a3 = -100. * this->zeta_3_A * (9. * this->omega_10_A + 8.) + 25200. * this->zeta_3_V - 48. * this->a_1_par + 240. * this->a_2_par;
    double a4 = 525. * this->zeta_3_A * this->omega_10_A - 14700. * this->zeta_3_V - 180. * this->a_2_par;
    return -u * (a1 + u * (a2 + u * (a3 + u * a4))) / 4. + this->delta_t_p * (9. * u - 1.5) + this->delta_t_m * 6. * u + 3. * (this->delta_t_p + this->delta_t_m) * log(1 - u);
}

complex_t BaseQCDfCalculator::F_V(double v, bool bar) {
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    return .75 * (
        BV::h(v, this->m_c_pole, this->mu_b) * (this->C_bar[WCoef::C2] + this->C_bar[WCoef::C4] + this->C_bar[WCoef::C6] + l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.)) 
      + BV::h(v, this->m_b_pole, this->mu_b) * (this->C_bar[WCoef::C3] + this->C_bar[WCoef::C4] + this->C_bar[WCoef::C6]) 
      + BV::h(v, 0., this->mu_b) * (this->C_bar[WCoef::C3] + 3. * this->C_bar[WCoef::C4] + 3. * this->C_bar[WCoef::C6] - l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.)) 
      - 8. / 27. * (this->C_bar[WCoef::C3] - this->C_bar[WCoef::C5] - 15. * this->C_bar[WCoef::C6])
    );
}

double BaseQCDfCalculator::L(double q2) {
    if (fpeq(q2, 0.0)) return 1.0;

    double mb2 = std::pow(this->m_b_PS, 2);
    return (q2 - mb2) * std::log(1. - q2 / mb2) / q2;
}

complex_t BaseQCDfCalculator::C_perp_0(double q2, double sign, bool bar) {
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);

    // printf("Y(q2 = %.3f) = %.4e + %.4e i\n", q2, real(Y(q2)), imag(Y(q2)));
    // printf("Y_u(q2 = %.3f) = %.4e + %.4e i\n", q2, real(Y_u(q2)), imag(Y_u(q2)));

    if (ff_tp == B_FF_Type::FULL) {
        if (fpeq(q2, 0.0)) return 0.0;
        return q2 * l_u * Y_u(q2) / (2. * this->m_b_PS * this->m_B);
    } else {
        if (fpeq(q2, 0.0)) return C7;
        return C7 + q2 * (Y(q2) + l_u * Y_u(q2)) / (2. * this->m_b_PS * this->m_B);
    }
}

complex_t BaseQCDfCalculator::C_par_0(double q2, double sign, bool bar) {
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);

    if (ff_tp == B_FF_Type::FULL) {
        return -this->m_B * l_u * Y_u(q2) / (2. * this->m_b_PS);
    } else {
        return -C7 - this->m_B * (Y(q2) + l_u * Y_u(q2)) / (2. * this->m_b_PS);
    }   
}

complex_t BaseQCDfCalculator::C_perp_f(double q2, double sign, bool bar) {
    if (ff_tp == B_FF_Type::FULL) return 0.0;

    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);

    // printf("C7 = %.4e + %.4e i\n", real(C7), imag(C7));
    // printf("L_b(q2 = %.3f) = %.4e\n", q2, this->L_b);
    // printf("L(q2 = %.3f) = %.4e\n", q2, L(q2));
    // printf("Delta_M = %.4e\n", this->Delta_M);

    return C7 * (-2. * this->L_b - L(q2) + this->Delta_M);
}

complex_t BaseQCDfCalculator::C_par_f(double q2, double sign, bool bar) {
    if (ff_tp == B_FF_Type::FULL) return 0.0;

    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 * (2. * std::log(this->m_b_PS / this->mu_b) + 2. * L(q2) + this->Delta_M);
}

complex_t BaseQCDfCalculator::C_perp_nf(double q2, bool bar) {
    double s_hat = q2 / (this->m_b_PS * this->m_b_PS);
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t F_27 = BV::f_27(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_27_u(s_hat, this->L_b) * l_u;

    if (fpeq(q2, 0.0))
        return -(this->C_bar[WCoef::C2] * F_27 + this->C[WCoef::C8] * BV::f_87(s_hat, this->L_b)) / iobs_qcdp->get_constants()->C_F;

    complex_t F_19 = BV::f_19_PS(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_19_u(s_hat, this->L_b) * l_u;
    complex_t F_29 = BV::f_29_PS(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_29_u(s_hat, this->L_b) * l_u;

    return -(
        this->C_bar[WCoef::C2] * F_27 
      + this->C[WCoef::C8] * BV::f_87(s_hat, this->L_b)
      + q2 / (2. * this->m_b_PS * this->m_B) * (
            (this->C_bar[WCoef::C2] + this->C_bar[WCoef::C1] / 3.) * F_29
          + 2. * this->C_bar[WCoef::C1] * F_19
          + this->C[WCoef::C8] * BV::f_89(s_hat)
        )
    ) / iobs_qcdp->get_constants()->C_F;
}

complex_t BaseQCDfCalculator::C_par_nf(double q2, bool bar) {
    double s_hat = q2 / (this->m_b_PS * this->m_b_PS);
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t F_27 = BV::f_27(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_27_u(s_hat, this->L_b) * l_u;
    complex_t F_19 = BV::f_19_PS(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_19_u(s_hat, this->L_b) * l_u;
    complex_t F_29 = BV::f_29_PS(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_29_u(s_hat, this->L_b) * l_u;

    return (
        this->C_bar[WCoef::C2] * F_27
      + this->C[WCoef::C8] * BV::f_87(s_hat, this->L_b)
      + this->m_B / (2 * this->m_b_PS) * (
            (this->C_bar[WCoef::C2] + this->C_bar[WCoef::C1] / 3.) * F_29
           + 2. * this->C_bar[WCoef::C1] * F_19
           + this->C[WCoef::C8] * BV::f_89(s_hat))
    ) / iobs_qcdp->get_constants()->C_F;
}

complex_t BaseQCDfCalculator::T_par_p_p_f(double u, double q2, bool bar) {
    return 2. * T_perp_p_p_f(u, q2, bar);
}

complex_t BaseQCDfCalculator::T_par_p_m_f(double u, double q2, bool bar) {
    return 2. * T_perp_p_m_f(u, q2, bar);
}

complex_t BaseQCDfCalculator::T_perp_p_p_f(double u, double q2, bool bar) {
    complex_t C7 = this->C[WCoef::C7] + this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return 2. * this->m_B / (1. - u) / this->E(q2) * C7;
}

complex_t BaseQCDfCalculator::T_perp_p_m_f(double u, double q2, bool bar) {
    complex_t C7 = this->C[WCoef::C7] - this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return 2. * this->m_B / (1. - u) / this->E(q2) * C7;
} 

complex_t BaseQCDfCalculator::T_perp_p_nf(double u, double q2, bool bar) {
    double E = this->E(q2);
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t t_perp_mc = t_perp(u, this->m_c_pole, q2, E);
    complex_t t_perp_mb = t_perp(u, this->m_b_PS, q2, E);
    complex_t t_perp_0 = t_perp(u, 0.0, q2, E);

    complex_t c8_term = /* fpeq(q2, 0.0) ? 0.0 : */ -4.0 * e_d * this->C[WCoef::C8] / (u + (1 - u) * q2 / (this->m_B * this->m_B));
    return c8_term + this->m_B / (2 * this->m_b_PS) * (
            e_u * (
                t_perp_mc * (this->C_bar[WCoef::C2] + this->C_bar[WCoef::C4] - this->C_bar[WCoef::C6] + l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.))
                - t_perp_0 * l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.)
            )
            + e_d * (
                t_perp_mb * (this->C_bar[WCoef::C3] + this->C_bar[WCoef::C4] - this->C_bar[WCoef::C6] - 4. * this->m_b_PS / this->m_B * this->C_bar[WCoef::C5])
                + t_perp_0 * this->C_bar[WCoef::C3]
            )
        );
}

complex_t BaseQCDfCalculator::T_par_p_nf(double u, double q2, bool bar) {
    double E = this->E(q2);
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t t_par_mc = t_par(u, this->m_c_pole, q2, E);
    complex_t t_par_mb = t_par(u, this->m_b_PS, q2, E);
    complex_t t_par_0 = t_par(u, 0., q2, E);

    // printf("t_par_mc = %.4e + %.4e i\n", real(t_par_mc), imag(t_par_mc));
    // printf("t_par_mb = %.4e + %.4e i\n", real(t_par_mb), imag(t_par_mb));
    // printf("t_par_0 = %.4e + %.4e i\n", real(t_par_0), imag(t_par_0));

    return this->m_B / this->m_b_PS * (
        e_u * (
            t_par_mc * (this->C_bar[WCoef::C2] + this->C_bar[WCoef::C4] - this->C_bar[WCoef::C6] + l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.))
            - t_par_0 * l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.)
        )
        + e_d * (
            t_par_mb * (this->C_bar[WCoef::C3] + this->C_bar[WCoef::C4] - this->C_bar[WCoef::C6])
            + t_par_0 * this->C_bar[WCoef::C3]
        )
    );
}

complex_t BaseQCDfCalculator::T_par_m_nf(double u, double q2, bool bar) {
    double v = this->m_B * this->m_B * (1 - u) + q2 * u;

    // printf("T_par_m C8 part = %.4e + %.4e i\n", real(8. * this->m_B * this->m_B * this->C[WCoef::C8] / v), imag(8. * this->m_B * this->m_B * this->C[WCoef::C8] / v));
    // printf("T_par_m F_V part = %.4e + %.4e i\n", real(8. * this->m_B / this->m_b_PS * F_V(v, bar)), imag(8. * this->m_B / this->m_b_PS * F_V(v, bar)));

    return 8. * this->m_B * this->m_B * this->C[WCoef::C8] / v + 8. * this->m_B / this->m_b_PS * F_V(v, bar);
}

complex_t BaseQCDfCalculator::inv_lambda_B_m(double q2) {
    double omega_0 = 2. * (this->m_B - this->m_b_PS) / 3.;
    double x = q2 / (this->m_B * omega_0);
    return std::exp(-x) / omega_0 * (I * PI - Ei(x));
}

complex_t BaseQCDfCalculator::I_perp_p(double q2, bool bar) {
    double pref = this->loop_f_mu_f / this->lambda_B_p;

    // printf("pref = %.4e\n", pref);
    // printf("T_perp_p_p_f(u = 0.5) = %.4e + %.4e i\n", real(T_perp_p_p_f(0.5, q2, bar)), imag(T_perp_p_p_f(0.5, q2, bar)));
    // printf("T_perp_p_nf(u = 0.5) = %.4e + %.4e i\n", real(T_perp_p_nf(0.5, q2, bar)), imag(T_perp_p_nf(0.5, q2, bar)));

    if (this->ff_tp == B_FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            return phi_X(u, this->a_1_perp, this->a_2_perp) * (T_perp_p_p_f(u, q2, bar) + T_perp_p_nf(u, q2, bar));
        };
        return pref * c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            return phi_X(u, this->a_1_perp, this->a_2_perp) * T_perp_p_nf(u, q2, bar);
        };
        return pref * c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BaseQCDfCalculator::I_perp_m(double q2, bool bar) {
    double pref = this->loop_f_mu_f / this->lambda_B_p;

    if (this->ff_tp == B_FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            return phi_X(u, this->a_1_perp, this->a_2_perp) * (T_perp_p_m_f(u, q2, bar) + T_perp_p_nf(u, q2, bar));
        };
        return pref * c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            return phi_X(u, this->a_1_perp, this->a_2_perp) * T_perp_p_nf(u, q2, bar);
        };
        return pref * c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BaseQCDfCalculator::I_par_p(double q2, bool bar) {
    if (this->ff_tp == B_FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * (T_par_p_p_f(u, q2, bar) + T_par_p_nf(u, q2, bar));
            complex_t i2 = phi * (T_par_m_0(bar) + this->loop_f_mu_f * T_par_m_nf(u, q2, bar));
            return this->loop_f_mu_f / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * T_par_p_nf(u, q2, bar);
            complex_t i2 = phi * (T_par_m_0(bar) + this->loop_f_mu_f * T_par_m_nf(u, q2, bar));
            return this->loop_f_mu_f / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BaseQCDfCalculator::I_par_m(double q2, bool bar) {
    if (this->ff_tp == B_FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * (T_par_p_m_f(u, q2, bar) + T_par_p_nf(u, q2, bar));
            complex_t i2 = phi * (T_par_m_0(bar) + this->loop_f_mu_f * T_par_m_nf(u, q2, bar));
            return this->loop_f_mu_f / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * T_par_p_nf(u, q2, bar);
            complex_t i2 = phi * (T_par_m_0(bar) + this->loop_f_mu_f * T_par_m_nf(u, q2, bar));
            return this->loop_f_mu_f / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_full, 0, 1, 1e-2);
    }
}