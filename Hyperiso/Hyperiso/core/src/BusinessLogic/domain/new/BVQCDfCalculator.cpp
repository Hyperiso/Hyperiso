#include "BVQCDfCalculator.h"
#include <mutex>

#include <chrono>
#include <iomanip>
#include <iostream>

namespace {
using qcdf_clock = std::chrono::steady_clock;

struct QcdfStat {
    long long calls {0};
    long long hits {0};
    double total_s {0.0};
    double max_s {0.0};
    mutable std::mutex mtx;

    void add(double dt, bool hit = false) {
        std::lock_guard<std::mutex> lock(mtx);
        ++calls;
        if (hit) ++hits;
        total_s += dt;
        if (dt > max_s) max_s = dt;
    }
};

struct QcdfScope {
    QcdfStat& stat;
    qcdf_clock::time_point t0 {qcdf_clock::now()};
    bool hit {false};

    explicit QcdfScope(QcdfStat& s) : stat(s) {}
    void mark_hit() { hit = true; }
    ~QcdfScope() {
        const double dt = std::chrono::duration<double>(qcdf_clock::now() - t0).count();
        stat.add(dt, hit);
    }
};

struct QcdfProfiler {
    QcdfStat ctor;
    QcdfStat precompute_static;
    QcdfStat T_triplet;
    QcdfStat T_perp_p;
    QcdfStat T_perp_m;
    QcdfStat T_par_m;
    QcdfStat triplet_C_perp_nf;
    QcdfStat triplet_C_perp_pm;
    QcdfStat triplet_xi_perp;
    QcdfStat triplet_I_perp_pm;
    QcdfStat triplet_C_par_nf;
    QcdfStat triplet_C_par_m;
    QcdfStat triplet_xi_par;
    QcdfStat triplet_E;
    QcdfStat triplet_I_par_m;
    QcdfStat delta_T_perp_WA;
    QcdfStat I_HSA_1;
    QcdfStat I_HSA_2;
    QcdfStat delta_T_perp_HSA;

    ~QcdfProfiler() {
        std::cerr << std::fixed << std::setprecision(6);
        std::cerr << "\n[BKPROF] BVQCDfCalculator cumulative timings\n";
        auto dump = [] (const char* name, const QcdfStat& s) {
            std::cerr << "[BKPROF] " << std::setw(24) << std::left << name
                      << " calls=" << std::setw(8) << s.calls
                      << " hits=" << std::setw(8) << s.hits
                      << " total=" << std::setw(10) << s.total_s << " s"
                      << " avg=" << (s.calls ? s.total_s / s.calls : 0.0) << " s"
                      << " max=" << s.max_s << " s\n";
        };
        dump("qcdf_ctor", ctor);
        dump("precompute_static", precompute_static);
        dump("T_triplet", T_triplet);
        dump("T_perp_p", T_perp_p);
        dump("T_perp_m", T_perp_m);
        dump("T_par_m", T_par_m);
        dump("triplet_C_perp_nf", triplet_C_perp_nf);
        dump("triplet_C_perp_pm", triplet_C_perp_pm);
        dump("triplet_xi_perp", triplet_xi_perp);
        dump("triplet_I_perp_pm", triplet_I_perp_pm);
        dump("triplet_C_par_nf", triplet_C_par_nf);
        dump("triplet_C_par_m", triplet_C_par_m);
        dump("triplet_xi_par", triplet_xi_par);
        dump("triplet_E", triplet_E);
        dump("triplet_I_par_m", triplet_I_par_m);
        dump("delta_T_perp_WA", delta_T_perp_WA);
        dump("I_HSA_1", I_HSA_1);
        dump("I_HSA_2", I_HSA_2);
        dump("delta_T_perp_HSA", delta_T_perp_HSA);
    }
};

QcdfProfiler& qcdf_prof() {
    static QcdfProfiler p;
    return p;
}
} // namespace

namespace {
inline complex_t get_or_zero(const std::map<WCoef, complex_t>& C, WCoef w) {
    auto it = C.find(w);
    return it == C.end() ? complex_t{} : it->second;
}
} // namespace

void BVQCDfCalculator::build_indexed_wilsons(const std::map<WCoef, complex_t>& C) {
    this->wc.C1 = get_or_zero(C, WCoef::C1);
    this->wc.C2 = get_or_zero(C, WCoef::C2);
    this->wc.C3 = get_or_zero(C, WCoef::C3);
    this->wc.C4 = get_or_zero(C, WCoef::C4);
    this->wc.C5 = get_or_zero(C, WCoef::C5);
    this->wc.C6 = get_or_zero(C, WCoef::C6);
    this->wc.C7 = get_or_zero(C, WCoef::C7);
    this->wc.C8 = get_or_zero(C, WCoef::C8);
    this->wc.C9 = get_or_zero(C, WCoef::C9);
    this->wc.C10 = get_or_zero(C, WCoef::C10);
    this->wc.CP7 = get_or_zero(C, WCoef::CP7);
    this->wc.CP9 = get_or_zero(C, WCoef::CP9);
    this->wc.CP10 = get_or_zero(C, WCoef::CP10);
    this->wc.CQ1 = get_or_zero(C, WCoef::CQ1);
    this->wc.CQ2 = get_or_zero(C, WCoef::CQ2);
    this->wc.CPQ1 = get_or_zero(C, WCoef::CPQ1);
    this->wc.CPQ2 = get_or_zero(C, WCoef::CPQ2);
}

BVQCDfCalculator::BVQCDfCalculator(int B_id, int V_id, double mu_b, const std::map<WCoef, complex_t> &C, std::shared_ptr<BVFFCalculator> ff_calculator, B_FF_Type ff_tp,
                        std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p,
                        std::shared_ptr<IObsQCDProxy> iobs_qcdp) :
    BaseQCDfCalculator(B_id, V_id, mu_b, C, ff_tp, p, iobs_qcdp)
{
    QcdfScope scope(qcdf_prof().ctor);
    this->ff_calculator = ff_calculator;
    build_indexed_wilsons(C);
    precompute_static_integrals();
}

void BVQCDfCalculator::precompute_static_integrals() {
    QcdfScope scope(qcdf_prof().precompute_static);
    if (this->static_integrals_ready) {
        scope.mark_hit();
        return;
    }

    auto iG_perp = [this] (double x) {
        double xbar = 1 - x;
        return phi_X(x, this->a_1_perp, this->a_2_perp) * BV::G(xbar, this->z_c) / (3 * xbar);
    };

    auto iH_perp = [this] (double x) {
        return gv_dga_4(x) * BV::G(1 - x, this->z_c);
    };

    auto iH2_perp = [this] (double x) {
        return BV::hard_kernel(1 - x, this->z_c) * phi_X(x, this->a_1_perp, this->a_2_perp);
    };

    double Nc = iobs_qcdp->get_constants()->Nc;
    const double T1_0 = ff_calculator->get(BV_FF::T1, 0.0);

    this->G_perp_cache = c_integrate(iG_perp, 0, 1, 1e-3);
    this->H_perp_cache = c_integrate(iH_perp, 0, 1, 1e-3);
    this->H_2_cache = -2 * PI2 * f_B * f_X_perp / (3 * Nc * m_B * lambda_B_p * T1_0) * c_integrate(iH2_perp, 0, 1, 1e-4);
    this->H_8_cache = 4. * PI2 * f_B * f_X_perp / (Nc * m_B * lambda_B_p * T1_0) * (1 - a_1_perp + a_2_perp);
    this->static_integrals_ready = true;
}


double BVQCDfCalculator::F_perp(double s) {
    if (fpeq(s, 0.0)) 
        return 1.0 + this->a_1_perp + this->a_2_perp;
    
    double d = s - 1.;
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double s2 = s * s;
    double ls = std::log(s);
    double f0 = (s + 1.) / d2 - 2. * s * ls / d3;
    double f1 = -(s2 + 10. * s + 1.) / d3 + 6. * s * (s + 1) * ls / d4;
    double f2 = (s + 1.) * (s2 + 28. * s + 1.) / d4 - 12. * s * (s2 + 3. * s + 1.) * ls / d5;
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
    double s3 = s2 * s;
    double s4 = s3 * s;
    double ls = std::log(s);
    double f0 = (s2 - 4 * s + 3 + 2. * ls) / d3;
    double f1 = -(s3 - 9 * s2 - 9. * s + 17. + 6. * (3. * s + 1.) * ls) / d4;
    double f2 = (-s4 + 16 * s3 + 108. * s2 - 80. * s - 43. - 12. * (6. * s2 + 8. * s + 1.) * ls) / d5;
    return f0 + this->a_1_perp * f1 + this->a_2_perp * f2;
}

complex_t BVQCDfCalculator::G_perp() {
    precompute_static_integrals();
    return this->G_perp_cache;
}

complex_t BVQCDfCalculator::H_perp() {
    precompute_static_integrals();
    return this->H_perp_cache;
}

complex_t BVQCDfCalculator::H_2() {
    precompute_static_integrals();
    return this->H_2_cache;
}

double BVQCDfCalculator::H_8() {
    precompute_static_integrals();
    return this->H_8_cache;
}


double BVQCDfCalculator::xi_perp_p_at(double q2) {
    if (fpeq(q2, 0.0)) {
        return this->ff_calculator->get(BV_FF::T1, 0.0) * (1 - this->loop_f_mu_b)
             + this->loop_f_mu_f * 12 * PI2 * this->f_B * this->f_X_perp
               / (iobs_qcdp->get_constants()->Nc * this->m_B * this->lambda_B_p) * F_perp(0.0);
    }
    return this->ff_calculator->get(BV_FF::XI_PERP, q2);
}

double BVQCDfCalculator::xi_perp_m_at(double q2) {
    if (fpeq(q2, 0.0)) {
        return this->ff_calculator->get(BV_FF::T1, 0.0) * (1 - this->alpha_s_mu_b * iobs_qcdp->get_constants()->C_F / (4 * PI))
             + this->alpha_s_mu_f * iobs_qcdp->get_constants()->C_F / (4 * PI)
               * 4 * PI2 * this->f_B * this->f_X_perp
               / (iobs_qcdp->get_constants()->Nc * this->m_B * this->lambda_B_p) * 3 * F_perp(0.0);
    }
    return this->ff_calculator->get(BV_FF::XI_PERP, q2);
}

BVQCDfCalculator::TTriplet BVQCDfCalculator::evaluate_T_triplet(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().T_triplet);

    TTriplet out {};

    complex_t c_perp_nf {};
    {
        QcdfScope s(qcdf_prof().triplet_C_perp_nf);
        c_perp_nf = C_perp_nf(q2, bar);
    }

    complex_t delta_perp_WA {};
    {
        QcdfScope s(qcdf_prof().delta_T_perp_WA);
        delta_perp_WA = delta_T_perp_WA(q2, bar);
    }
    const complex_t delta_perp_HSA = delta_T_perp_HSA(q2, bar);
    const complex_t delta_perp = delta_perp_WA + delta_perp_HSA;

    complex_t c_perp_p {};
    complex_t c_perp_m {};
    {
        QcdfScope s(qcdf_prof().triplet_C_perp_pm);
        c_perp_p = C_perp_0(q2, 1, bar) + this->loop_f_mu_b * (C_perp_f(q2, 1, bar) + c_perp_nf);
        c_perp_m = C_perp_0(q2, -1, bar) + this->loop_f_mu_b * (C_perp_f(q2, -1, bar) + c_perp_nf);
    }

    double xi_perp_p {};
    double xi_perp_m {};
    {
        QcdfScope s(qcdf_prof().triplet_xi_perp);
        xi_perp_p = xi_perp_p_at(q2);
        xi_perp_m = xi_perp_m_at(q2);
    }

    complex_t i_perp_p {};
    complex_t i_perp_m {};
    {
        QcdfScope s(qcdf_prof().triplet_I_perp_pm);
        i_perp_p = I_perp_p(q2, bar);
        i_perp_m = I_perp_m(q2, bar);
    }

    out.perp_p = xi_perp_p * c_perp_p + this->pref_perp * i_perp_p + delta_perp;
    out.perp_m = xi_perp_m * c_perp_m + this->pref_perp * i_perp_m + delta_perp;

    complex_t c_par_nf {};
    {
        QcdfScope s(qcdf_prof().triplet_C_par_nf);
        c_par_nf = C_par_nf(q2, bar);
    }

    complex_t c_par_m {};
    {
        QcdfScope s(qcdf_prof().triplet_C_par_m);
        c_par_m = C_par_0(q2, -1, bar) + this->loop_f_mu_b * (C_par_f(q2, -1, bar) + c_par_nf);
    }

    double xi_par {};
    {
        QcdfScope s(qcdf_prof().triplet_xi_par);
        xi_par = this->ff_calculator->get(BV_FF::XI_PAR, q2);
    }

    double e_q2 {};
    {
        QcdfScope s(qcdf_prof().triplet_E);
        e_q2 = this->E(q2);
    }

    complex_t i_par_m {};
    {
        QcdfScope s(qcdf_prof().triplet_I_par_m);
        i_par_m = I_par_m(q2, bar);
    }

    out.par_m = xi_par * c_par_m + this->pref_par * this->m_X / e_q2 * i_par_m;

    return out;
}

complex_t BVQCDfCalculator::T_perp_p(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().T_perp_p);
    complex_t C_perp_p = C_perp_0(q2, 1, bar) + this->loop_f_mu_b * (C_perp_f(q2, 1, bar) + C_perp_nf(q2, bar));
    return xi_perp_p_at(q2) * C_perp_p + this->pref_perp * I_perp_p(q2, bar) + delta_T_perp_WA(q2, bar) + delta_T_perp_HSA(q2, bar);
}

complex_t BVQCDfCalculator::T_perp_m(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().T_perp_m);
    complex_t C_perp_m = C_perp_0(q2, -1, bar) + this->loop_f_mu_b * (C_perp_f(q2, -1, bar) + C_perp_nf(q2, bar));
    return xi_perp_m_at(q2) * C_perp_m + this->pref_perp * I_perp_m(q2, bar) + delta_T_perp_WA(q2, bar) + delta_T_perp_HSA(q2, bar);
}

complex_t BVQCDfCalculator::T_par_p(double q2, bool bar) {
    complex_t C_par_p = C_par_0(q2, 1, bar) + this->loop_f_mu_b * (C_par_f(q2, 1, bar) + C_par_nf(q2, bar));
    return this->ff_calculator->get(BV_FF::XI_PAR, q2) * C_par_p + this->pref_par * this->m_X / this->E(q2) * I_par_p(q2, bar);
}

complex_t BVQCDfCalculator::T_par_m(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().T_par_m);
    complex_t C_par_m = C_par_0(q2, -1, bar) + this->loop_f_mu_b * (C_par_f(q2, -1, bar) + C_par_nf(q2, bar));
    return this->ff_calculator->get(BV_FF::XI_PAR, q2) * C_par_m + this->pref_par * this->m_X / this->E(q2) * I_par_m(q2, bar);
}

complex_t BVQCDfCalculator::Delta_par(double q2) {
    return 1. + 2 * this->loop_f_mu_b * (
        L(q2) - 1 - 3. * PI2 * q2 * this->f_B * this->f_X_par * this->m_X / (iobs_qcdp->get_constants()->Nc * this->m_B * this->lambda_B_p * this->ff_calculator->get(BV_FF::XI_PAR, q2) * std::pow(this->E(q2), 3)) * F_perp(0.0)
    );
}

complex_t BVQCDfCalculator::I_HSA_1(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().I_HSA_1);
    auto key = std::make_pair(bar, q2);
    auto it = this->i_hsa1_cache.find(key);
    if (it != this->i_hsa1_cache.end()) {
        scope.mark_hit();
        return it->second;
    }

    auto f = [q2, bar, this] (double u) {
        double phi_u = phi_X(u, this->a_1_perp, this->a_2_perp);
        double v = this->m_B * this->m_B * (1 - u) + u * q2;
        return phi_u * this->m_B * this->m_B / v * F_V(v, bar);
    };

    complex_t value = c_integrate(f, 0, 1, 1e-3);
    this->i_hsa1_cache.emplace(key, value);
    return value;
}

complex_t BVQCDfCalculator::I_HSA_2(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().I_HSA_2);
    auto key = std::make_pair(bar, q2);
    auto it = this->i_hsa2_cache.find(key);
    if (it != this->i_hsa2_cache.end()) {
        scope.mark_hit();
        return it->second;
    }

    auto f = [q2, bar, this] (double u) {
        double int_phi_par = gv_dga_4(u);
        double v = this->m_B * this->m_B * (1 - u) + u * q2;
        return int_phi_par * F_V(v, bar);
    };

    complex_t value = c_integrate(f, 0, 1, 1e-2);
    this->i_hsa2_cache.emplace(key, value);
    return value;
}

complex_t BVQCDfCalculator::delta_T_perp_WA(double q2, bool bar) {
    double pref = this->e_q * 2. * PI2 * this->f_B / (this->m_b_PS * this->m_B);
    complex_t l_u = bar ? std::conj(this->lambda_hat_u) : this->lambda_hat_u;
    complex_t W_perp = this->wc.C3 + 4. / 3. * (this->wc.C4 + 3. * this->wc.C5 + 4. * this->wc.C6);
    complex_t W_par = this->wc.C3 + 4. / 3. * this->wc.C4 + 16. * this->wc.C5 + 64. / 3. * this->wc.C6;
    W_par += delta_qu * l_u * -3. * this->wc.C2;
    if (this->src_block == "B_phi") {
        W_par += -l_u * (4. / 3 * this->wc.C1 + this->wc.C2);
        W_par += 12. * (this->wc.C3 + 10. * this->wc.C5);
    }
    
    double s_hat = q2 / (this->m_B * this->m_B);

    // printf("pref_1(s = %.3f) = %.4e\n", q2, pref * -2. * this->f_X_perp);
    // printf("W_perp(s = %.3f) = %.4e + %.4e i\n", q2, real(W_perp), imag(W_perp));
    // printf("F_perp(s = %.3f) = %.4e\n", q2, F_perp(s_hat));

    return pref * (
        -2. * this->f_X_perp * W_perp * F_perp(s_hat)
      + this->f_X_par * this->m_X * W_par / (3. * (1 - s_hat) * this->lambda_B_p)
    );
}

complex_t BVQCDfCalculator::delta_T_perp_HSA(double q2, bool bar) {
    QcdfScope scope(qcdf_prof().delta_T_perp_HSA);
    double pref = this->e_q * this->loop_f_mu_f * 4 * PI2 * this->f_B / (iobs_qcdp->get_constants()->Nc * this->m_b_PS * this->m_B);
    double s_hat = q2 / (this->m_B * this->m_B);

    // printf("X_perp = %.4e + %.4e i\n", std::real(X_perp(s_hat)), std::imag(X_perp(s_hat)));
    // printf("I_HSA_1 = %.4e + %.4e i\n", std::real(I_HSA_1(q2, bar)), std::imag(I_HSA_1(q2, bar)));
	// printf("I_HSA_2 = %.4e + %.4e i\n", std::real(I_HSA_2(q2, bar)), std::imag(I_HSA_2(q2, bar)));
    // printf("pref_C8_Xperp = %.4e\n", pref * 3. * this->m_b_PS / this->m_B * this->f_X_perp);
    // printf("pref_I1 = %.4e\n", pref * 2. * this->f_X_perp);
    // printf("pref_I2 = %.4e\n", -pref * this->m_X * this->f_X_par / ((1 - s_hat) * this->lambda_B_p));

    return pref * (
        3. * this->wc.C8 * this->m_b_PS / this->m_B * this->f_X_perp * X_perp(s_hat)
      + 2. * this->f_X_perp * I_HSA_1(q2, bar)
      - this->m_X * this->f_X_par / ((1 - s_hat) * this->lambda_B_p) * I_HSA_2(q2, bar)
    );
}