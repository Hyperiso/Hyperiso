#include "BVQCDfCalculator.h"

BVQCDfCalculator::BVQCDfCalculator(int B_id, int V_id, double mu_b, const std::map<WCoef, complex_t> &C, std::shared_ptr<BVFFCalculator> ff_calculator, B_FF_Type ff_tp) : 
    BaseQCDfCalculator(B_id, V_id, mu_b, C, ff_tp)
{
    this->ff_calculator = ff_calculator;
}

double BVQCDfCalculator::F_perp(double s) {
    if (fpeq(s, 0.0)) 
        return 1.0 + this->a_1_perp + 2.0 * this->a_2_perp;
    
    double d = s - 1.;
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double s2 = s * s;
    double ls = std::log(s);
    double f0 = (s + 1.) / d2 - 2. * s * ls / d3;
    double f1 = -(s2 + 10. * s + 1.) / d3 + 6. * s * (s + 1) * ls / d4;
    double f2 = (s + 1.) + (s2 + 28. * s + 1.) / d4 - 12. * s * (s2 + 3. * s + 1.) * ls / d5;
    return f0 + this->a_1_perp * f1 + this->a_2_perp * f2;
}

double BVQCDfCalculator::X_perp(double s) {
    if (fpeq(s, 0.0)) {
        double cutoff = this->Lambda_h / this->m_B;
        return -2 * (1 + 3 * this->a_1_perp + 6 * this->a_2_perp) * log(cutoff) - (1 + 11 * this->a_1_perp + 31 * this->a_2_perp) + 12 * cutoff * (this->a_1_perp + 5 * this->a_2_perp);
    }

    double d = s - 1;
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double s2 = s * s;
    double ls = std::log(s);
    double f0 = (s - 3.) / d2 + 2. * ls / d3;
    double f1 = -(s2 - 8. * s - 17.) / d3 + 6. * (3. * s + 1.) * ls / d4;
    double f2 = -(s * s2 - 15. * s2 - 123. * s - 43.) / d4 - 12. * (6. * s2 + 8. * s + 1.) * ls / d5;
    return f0 + this->a_1_perp * f1 + this->a_2_perp * f2;
}

complex_t BVQCDfCalculator::G_perp() {
    auto iG_perp = [this] (double x) {
        double xbar = 1 - x;
        return phi_X(x, this->a_1_perp, this->a_2_perp) * BV::G(xbar, this->z_c) / (3 * xbar);
    };

    return c_integrate(iG_perp, 0, 1, 1e-3);
}

complex_t BVQCDfCalculator::H_perp() {
    auto iH_perp = [this] (double x) {
        return gv_dga_4(x) * BV::G(1 - x, this->z_c);
    }; 
    return c_integrate(iH_perp, 0, 1, 1e-3);
}

complex_t BVQCDfCalculator::H_2() {
    auto iH2_perp = [this] (double x) {
        return BV::hard_kernel(1 - x, this->z_c) * phi_X(x, this->a_1_perp, this->a_2_perp);
    }; 

    double Nc = ObsQCDProxy().get_constants()->Nc;
    return -2 * PI2 * f_B * f_X_perp / (3 * Nc * m_B * lambda_B_p * ff_calculator->get(BV_FF::T1, 0.0)) * c_integrate(iH2_perp, 0, 1, 1e-4);
}

double BVQCDfCalculator::H_8() {
    double Nc = ObsQCDProxy().get_constants()->Nc;
    return 4 * PI2 * f_B * f_X_perp / (Nc * m_B * lambda_B_p * ff_calculator->get(BV_FF::T1, 0.0)) * (1 - a_1_perp + a_2_perp);
}

complex_t BVQCDfCalculator::T_perp_p(double q2, bool bar) {
    complex_t C_perp_p;

    if (this->ff_tp == B_FF_Type::SOFT) {
        C_perp_p = C_perp_0(q2, 1, bar) + this->alpha_s_mu_b / (4. * PI) * (C_perp_f(q2, 1, bar) + C_perp_nf(q2, bar));
    } else {
        C_perp_p = this->alpha_s_mu_b / (4. * PI) * C_perp_nf(q2, bar);
    }

    return this->ff_calculator->get(BV_FF::XI_PERP, q2) * C_perp_p + this->pref_perp * I_perp_p(q2, bar) + delta_T_perp_WA(q2) + delta_T_perp_HSA(q2, bar);
}

complex_t BVQCDfCalculator::T_perp_m(double q2, bool bar) {
    complex_t C_perp_m;

    if (this->ff_tp == B_FF_Type::SOFT) {
        C_perp_m = C_perp_0(q2, -1, bar) + this->alpha_s_mu_b / (4. * PI) * (C_perp_f(q2, -1, bar) + C_perp_nf(q2, bar));
    } else {
        C_perp_m = this->alpha_s_mu_b / (4. * PI) * C_perp_nf(q2, bar);
    }

    return this->ff_calculator->get(BV_FF::XI_PERP, q2) * C_perp_m + this->pref_perp * I_perp_m(q2, bar) + delta_T_perp_WA(q2) + delta_T_perp_HSA(q2, bar);
}

complex_t BVQCDfCalculator::T_par_p(double q2, bool bar) {
    complex_t C_par_p;

    if (this->ff_tp == B_FF_Type::SOFT) {
        C_par_p = C_par_0(q2, 1, bar) + this->alpha_s_mu_b / (4. * PI) * (C_par_f(q2, 1, bar) + C_par_nf(q2, bar));
    } else {
        C_par_p = this->alpha_s_mu_b / (4. * PI) * C_par_nf(q2, bar);
    }

    return this->ff_calculator->get(BV_FF::XI_PAR, q2) * C_par_p + this->pref_par / this->E(q2) * I_par_p(q2, bar);
}

complex_t BVQCDfCalculator::T_par_m(double q2, bool bar) {
    complex_t C_par_m;

    if (this->ff_tp == B_FF_Type::SOFT) {
        C_par_m = C_par_0(q2, -1, bar) + this->alpha_s_mu_b / (4. * PI) * (C_par_f(q2, -1, bar) + C_par_nf(q2, bar));
    } else {
        C_par_m = this->alpha_s_mu_b / (4. * PI) * C_par_nf(q2, bar);
    }

    return this->ff_calculator->get(BV_FF::XI_PAR, q2) * C_par_m + this->pref_par / this->E(q2) * I_par_m(q2, bar);
}

complex_t BVQCDfCalculator::Delta_par(double q2) {
    return 1. + this->alpha_s_mu_b * ObsQCDProxy().get_constants()->C_F / (2. * PI) * (
        L(q2) - 1 - 3. * PI2 * q2 * this->f_B * this->f_X_par * this->m_X / (ObsQCDProxy().get_constants()->Nc * this->m_B * this->lambda_B_p * this->ff_calculator->get(BV_FF::XI_PAR, q2) * std::pow(this->E(q2), 3)) * F_perp(0.0)
    );
}

complex_t BVQCDfCalculator::I_HSA_1(double q2, bool bar) {
    auto f = [q2, bar, this] (double u) {
        double phi_u = phi_X(u, this->a_1_par, this->a_2_par);
        double v = this->m_B * this->m_B * (1 - u) + u * q2;
        return phi_u * this->m_B * this->m_B / v * F_V(v, bar);
    };
    
    return c_integrate(f, 0, 1, 1e-2);
}

complex_t BVQCDfCalculator::I_HSA_2(double q2, bool bar) {
    auto f = [q2, bar, this] (double u) {
        double int_phi_par = gv_dga_4(u);
        double v = this->m_B * this->m_B * (1 - u) + u * q2;
        return int_phi_par * F_V(v, bar);
    };
    
    return c_integrate(f, 0, 1, 1e-2);
}

complex_t BVQCDfCalculator::delta_T_perp_WA(double q2) {
    double pref = this->e_q * 2. * PI2 * this->f_B / (this->m_b_PS * this->m_B);
    complex_t W_perp = this->C[WCoef::C3] + 4. / .3 * (this->C[WCoef::C4] + 3. * this->C[WCoef::C5] + 4. * this->C[WCoef::C6]);
    complex_t W_par = this->C_bar[WCoef::C3] + 3. * this->C_bar[WCoef::C4];
    W_par += delta_qu * -3. * this->C[WCoef::C2];
    double s_hat = q2 / (this->m_B * this->m_B);
    return pref * (
        -2. * this->f_X_perp * W_perp * F_perp(s_hat)
      + this->f_X_par * this->m_X * W_par / (3. * (1 - s_hat) * this->lambda_B_p)
    );
}

complex_t BVQCDfCalculator::delta_T_perp_HSA(double q2, bool bar) {
    double pref = this->e_q * this->alpha_s_mu_b * ObsQCDProxy().get_constants()->C_F * PI * this->f_B / (ObsQCDProxy().get_constants()->Nc * this->m_b_PS * this->m_B);
    double s_hat = q2 / (this->m_B * this->m_B);
    return pref * (
        3. * this->C[WCoef::C8] * this->m_b_PS / this->m_B * this->f_X_perp * X_perp(s_hat)
      + 2. * this->f_X_perp * I_HSA_1(q2, bar)
      - this->m_X * this->f_X_par / ((1 - s_hat) * this->lambda_B_p) * I_HSA_2(q2, bar)
    );
}