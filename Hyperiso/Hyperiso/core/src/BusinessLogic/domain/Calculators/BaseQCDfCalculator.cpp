#include "BaseQCDfCalculator.h"

double BaseQCDfCalculator::E(double q2) {
    return (m_B * m_B + m_X * m_X - q2) / (2 * m_B);
}

BaseQCDfCalculator::BaseQCDfCalculator(int B_id, int X_id, double mu_b, const std::map<WCoef, complex_t> &C, B_FF_Type ff_tp) :
    mu_b(mu_b), C(C), ff_tp(ff_tp)
{
    if (!this->allowed_decays.contains({B_id, X_id})) {
        LOG_ERROR("ValueError", "Wrong meson PDG code in BaseQCDfCalculator constructor:", B_id, ",", X_id);
    }

    this->src_block = allowed_decays.at({B_id, X_id});
    this->delta_qu = double (B_id == 521);
    this->fill_wilson_bar_cache();

    ObsParameterProxy p;
    double beta_0 = ObsQCDProxy().get_constants()->beta[5][0]; // TODO : Link with get_nf vs. hard-coded ?
    auto run = [this, beta_0] (double value_1gev, double eta, double gamma) { return value_1gev * pow(eta, gamma / beta_0); };

    this->Lambda_h = p(ParamId{ParameterType::DECAY, this->src_block, 14});
    double mu_f = sqrt(this->mu_b * this->Lambda_h);
    this->alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(this->mu_b, MassType::POLE, MassType::POLE));
    this->alpha_s_mu_f = ObsQCDProxy()(AlphasConfig(mu_f, MassType::POLE, MassType::POLE));
    this->m_c_pole = p(ParamId{ParameterType::SM, "QCD", 4});
    this->m_b_pole = p(ParamId{ParameterType::SM, "QCD", {5, 2}});
    double eta_f = this->alpha_s_mu_f / ObsQCDProxy()(AlphasConfig(1.0, MassType::POLE, MassType::POLE));
    this->m_b_PS = this->m_b_pole - 4 * ObsQCDProxy()(AlphasConfig(this->m_b_pole, MassType::POLE, MassType::POLE)) * mu_f / (3 * PI);
    this->m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", B_id});
    this->m_X = p(ParamId{ParameterType::FLAVOR, "FMASS", X_id});
    this->f_B = p(ParamId{ParameterType::FLAVOR, "FCONST", {B_id, 1}});
    this->f_X_par = p(ParamId{ParameterType::FLAVOR, "FCONST", {X_id, 1}});
    this->lambda_B_p = p(ParamId{ParameterType::DECAY, this->src_block, 13}) / (1. - this->alpha_s_mu_f * log(pow(mu_f, 2)) * 1.8 / (3. * PI));
    this->lambda_hat_u = std::conj(p(ParamId{ParameterType::SM, "VCKM", {0, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {0, 2}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}));
    this->a_1_par = run(p(ParamId{ParameterType::DECAY, this->src_block, {8, 1}}), eta_f, gamma_par(1));
    this->a_2_par = run(p(ParamId{ParameterType::DECAY, this->src_block, {8, 2}}), eta_f, gamma_par(2));
    this->e_q = B_id == 321 ? e_u : e_d;
    this->z_c = std::pow(this->m_c_pole / this->m_b_PS, 2);
    this->L_b = std::log(this->mu_b / this->m_b_PS);
    this->Delta_M = -6. * this->L_b - 4. * (1. - mu_f / this->m_b_PS);
    int Nc = ObsQCDProxy().get_constants()->Nc;
    this->pref_par = PI2 * this->f_B * this->f_B * this->f_X_par / (Nc * this->m_X);
    this->T_par_m_0 = 4. * this->m_B / this->m_b_PS * (3. * this->lambda_hat_u * this->delta_qu * this->C[WCoef::C2] - this->C_bar[WCoef::C3] - 3. * this->C_bar[WCoef::C4]);

    bool isV = (X_id == 313 || X_id == 323 || X_id == 331);
    if (isV) {
        this->f_X_perp = run(p(ParamId{ParameterType::FLAVOR, "FCONST", {X_id, 2}}), eta_f, ObsQCDProxy().get_constants()->C_F);
        this->a_1_perp = run(p(ParamId{ParameterType::DECAY, this->src_block, {7, 1}}), eta_f, gamma_perp(1));
        this->a_2_perp = run(p(ParamId{ParameterType::DECAY, this->src_block, {7, 2}}), eta_f, gamma_perp(2));
        this->zeta_3_A = p(ParamId{ParameterType::DECAY, this->src_block, 9});
        this->zeta_3_V = p(ParamId{ParameterType::DECAY, this->src_block, 10});
        this->omega_10_A = p(ParamId{ParameterType::DECAY, this->src_block, 11});
        this->delta_t_p = p(ParamId{ParameterType::DECAY, this->src_block, {12, 1}});
        this->delta_t_m = p(ParamId{ParameterType::DECAY, this->src_block, {12, 2}});
        this->pref_perp = PI2 * this->f_B * this->f_X_perp / (Nc * this->m_B);
    }

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
    double mB2 = this->m_B * this->m_B;
    if(fpeq(q2, 0.)) {
		if (fpeq(m_q, 0.)) return 4./(1.-u);
		double epsilon=1.e-10;
		complex_t xp=0.5+std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*mB2));
		complex_t xm=0.5-std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*mB2));
		return 4./(1.-u)*(1.+2.*m_q*m_q/(1.-u)/mB2*(BV::L_1(xp)+BV::L_1(xm)));
	} else {
        double s_hat = q2 / mB2;
        double mq_hat = m_q / this->m_B;
        return 2.*this->m_B/(1.-u)/E_Kstar*BV::I_1(u,s_hat,mq_hat)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(BV::B_0((1.-u)*mB2+u*q2,m_q)-BV::B_0(q2,m_q));
    }
}

complex_t BaseQCDfCalculator::t_par(double u, double m_q, double q2, double E_Kstar) {
    double mB2 = this->m_B * this->m_B;
    double s_hat = q2 / mB2;
    double mq_hat = m_q / this->m_B;
    return 2.*this->m_B/(1.-u)/E_Kstar*BV::I_1(u,s_hat,mq_hat)+((1.-u)*mB2+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(BV::B_0((1.-u)*mB2+u*q2,m_q)-BV::B_0(q2,m_q));
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
    return (q2 - mb2) * std::log(1 - q2 / mb2) / q2;
}

complex_t BaseQCDfCalculator::C_perp_0(double q2, double sign, bool bar) {
    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);

    if (fpeq(q2, 0.0)) return C7;

    return C7 + q2 * (Y(q2) + this->lambda_hat_u * Y_u(q2)) / (2. * this->m_b_PS * this->m_B);
}

complex_t BaseQCDfCalculator::C_par_0(double q2, double sign, bool bar) {
    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 - this->m_B * (Y(q2) + this->lambda_hat_u * Y_u(q2)) / (2. * this->m_b_PS);
}

complex_t BaseQCDfCalculator::C_perp_f(double q2, double sign, bool bar) {
    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return C7 * (2. * std::log(this->m_b_PS / this->mu_b) - L(q2) + this->Delta_M);
}

complex_t BaseQCDfCalculator::C_par_f(double q2, double sign, bool bar) {
    complex_t C7 = this->C[WCoef::C7] + sign * this->C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 * (2. * std::log(this->m_b_PS / this->mu_b) + 2. * L(q2) + this->Delta_M);
}

complex_t BaseQCDfCalculator::C_perp_nf(double q2, bool bar) {
    double s_hat = q2 / (this->m_b_PS * this->m_b_PS);
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    BV::f_27(s_hat, this->L_b, this->z_c);
    BV::f_27_u(s_hat, this->L_b);
    complex_t F_27 = BV::f_27(s_hat, this->L_b, this->z_c) * (1. + l_u) + BV::f_27_u(s_hat, this->L_b) * l_u;

    if (fpeq(q2, 0.0))
        return -(this->C_bar[WCoef::C2] * F_27 + this->C[WCoef::C8] * BV::f_87(s_hat, this->L_b)) / ObsQCDProxy().get_constants()->C_F;

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
    ) / ObsQCDProxy().get_constants()->C_F;
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
    ) / ObsQCDProxy().get_constants()->C_F;
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
    complex_t c8_term = -4.0 * e_d * this->C[WCoef::C8] / (u + (1 - u) * q2 / (this->m_B * this->m_B));
    return c8_term + this->m_B / (2 * this->m_b_PS) * (
            e_u * (
                t_perp_mc * (this->C_bar[WCoef::C2] + this->C_bar[WCoef::C4] - this->C_bar[WCoef::C6] + l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.))
                - t_perp_0 * l_u * (this->C[WCoef::C2] - this->C[WCoef::C1] / 6.)
            )
            + e_d * (
                t_perp_mb * (this->C_bar[WCoef::C3] + this->C_bar[WCoef::C4] - this->C_bar[WCoef::C6] - 4 * this->m_b_PS / this->m_B * this->C_bar[WCoef::C5])
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
    return 8. * this->m_B * this->m_B * this->C[WCoef::C8] / v + 8. * this->m_B / this->m_b_PS * F_V(v, bar);
}

complex_t BaseQCDfCalculator::inv_lambda_B_m(double q2) {
    double omega_0 = 2. * (this->m_B - this->m_b_PS) / 3.;
    double x = q2 / (this->m_B * omega_0);
    return std::exp(-x) / omega_0 * (I * PI - Ei(x));
}

complex_t BaseQCDfCalculator::I_perp_p(double q2, bool bar) {
    double pref = this->alpha_s_mu_f / (4. * PI) * ObsQCDProxy().get_constants()->C_F / this->lambda_B_p;

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
    double pref = this->alpha_s_mu_f / (4. * PI) * ObsQCDProxy().get_constants()->C_F / this->lambda_B_p;

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
            double fact = this->alpha_s_mu_f * ObsQCDProxy().get_constants()->C_F / (4 * PI);
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * (T_par_p_p_f(u, q2, bar) + T_par_p_nf(u, q2, bar));
            complex_t i2 = phi * (this->T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            double fact = this->alpha_s_mu_f * ObsQCDProxy().get_constants()->C_F / (4 * PI);
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * T_par_p_nf(u, q2, bar);
            complex_t i2 = phi * (this->T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BaseQCDfCalculator::I_par_m(double q2, bool bar) {
    if (this->ff_tp == B_FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            double fact = this->alpha_s_mu_f * ObsQCDProxy().get_constants()->C_F / (4 * PI);
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * (T_par_p_m_f(u, q2, bar) + T_par_p_nf(u, q2, bar));
            complex_t i2 = phi * (this->T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            double fact = this->alpha_s_mu_f * ObsQCDProxy().get_constants()->C_F / (4 * PI);
            double phi = phi_X(u, this->a_1_par, this->a_2_par);
            complex_t i1 = phi * T_par_p_nf(u, q2, bar);
            complex_t i2 = phi * (this->T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / this->lambda_B_p * i1 + this->e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_full, 0, 1, 1e-2);
    }
}