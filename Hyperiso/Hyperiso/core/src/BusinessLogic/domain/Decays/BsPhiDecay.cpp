#include <algorithm>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

#include "BsPhiDecay.h"

void BsPhiDecay::load_params() {
    fill_wilson_cache();

    cache.ff_calculator = BVFFCalculator(531, 333, p, cfg.ff_src);
    cache.mu_b = (*p)(ParamId{ParameterType::WILSON, "B_SCALE", 1}, DataType::VALUE);

    cache.qcdf_calculator = BVQCDfCalculator(
        531, 333,
        cache.mu_b,
        cache.C,
        std::make_shared<BVFFCalculator>(cache.ff_calculator),
        cfg.ff_type,
        p,
        iobs_qcdp
    );

    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.G_F = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    cache.m_s = (*p)(ParamId{ParameterType::SM, "MASS", 3}, DataType::VALUE);
    cache.alpha_s_mu_b = (*iobs_qcdp)(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.m_c_m_c = (*p)(ParamId{ParameterType::SM, "MASS", 4}, DataType::VALUE);
    cache.m_b_mu_b = (*iobs_qcdp)(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_b_m_b = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 5}, DataType::VALUE);
    cache.m_b_PS = (*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE) - 4 * (*iobs_qcdp)(AlphasConfig((*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE), MassType::POLE, MassType::POLE)) * sqrt(cache.mu_b * (*p)(ParamId{ParameterType::DECAY, "B_phi", 14}, DataType::VALUE)) / (3 * PI);
    cache.L_b = std::log(cache.mu_b / cache.m_b_PS);
    cache.m_Bs = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 531}, DataType::VALUE);
    cache.m_phi = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 333}, DataType::VALUE);
    cache.life_Bs = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 531}, DataType::VALUE) / HBAR;
    cache.lambda_hat_u = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE) 
                            / (std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE));
    cache.kappa = 1 - 2. * cache.alpha_s_mu_b / (3. * PI) * std::log(cache.mu_b / cache.m_b_m_b);
    cache.Delta_M = -6. * cache.L_b - 4. * (1 - sqrt(cache.mu_b * (*p)(ParamId{ParameterType::DECAY, "B_phi", 14}, DataType::VALUE)) / cache.m_b_PS);
    cache.N_0 = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * cache.G_F * cache.alpha_em / (std::sqrt(3072. * std::pow(PI, 5) * std::pow(cache.m_Bs, 3)));
    cache.q2_max = std::pow(cache.m_Bs - cache.m_phi, 2);
    cache.q2_low = (*p)(ParamId{ParameterType::DECAY, "B_phi", {15, 1}}, DataType::VALUE);
    cache.q2_high = (*p)(ParamId{ParameterType::DECAY, "B_phi", {15, 2}}, DataType::VALUE);
    cache.ys = (*p)(ParamId{ParameterType::DECAY, "B_ll", 1}, DataType::VALUE);
    // ASK : hardcoded in SI
    // double beta_s = std::arg(-std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}})) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}) 
    //                         / (std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {1, 2}})) * (*p)(ParamId{ParameterType::SM, "VCKM", {1, 1}})));
    // cache.phi_s = 2 * beta_s;
    // cache.up = std::exp(I * cache.phi_s);
    // cache.um = std::exp(-I * cache.phi_s);

    cache.phi_s = 0.04;
    cache.up = 1. + cache.phi_s * I;
    cache.um = 1. - cache.phi_s * I;

    complex_t eipi4 = std::exp(I * PI / 4.0);

    for (size_t i = 0; i < 6; i++) {
        cache.A_had_err_low_0[i] = (*p)(ParamId{ParameterType::DECAY, "B_phi", {18, 1, i + 1}}, DataType::VALUE);
        cache.A_had_err_low_1[i] = (*p)(ParamId{ParameterType::DECAY, "B_phi", {18, 2, i + 1}}, DataType::VALUE);
    }

    for (size_t i = 0; i < 8; i++) {
        cache.A_had_err_high[i] = (*p)(ParamId{ParameterType::DECAY, "B_phi", {18, 3, i + 1}}, DataType::VALUE);
    }

    load_cfg_dep_params();

    // #########################################
    // #                DEBUG                  #
    // #########################################

    // printf("y_s = %.4e\n", cache.ys);
	// printf("up = %.4e + %.4e i\n", real(cache.up), imag(cache.up));
	// printf("um = %.4e + %.4e i\n", real(cache.um), imag(cache.um));

    // double q2 = 15.0; 
    // double u = 0.5;
    // complex_t Tperpp = cache.qcdf_calculator.T_perp_p(q2, false);
    // complex_t Tperpm = cache.qcdf_calculator.T_perp_m(q2, false);
    // complex_t Tparm = cache.qcdf_calculator.T_par_m(q2, false);

    // printf("T_perp_p(s = %.3f) = %.4e + %.4e i\n", q2, std::real(Tperpp), std::imag(Tperpp));
    // printf("T_perp_m(s = %.3f) = %.4e + %.4e i\n", q2, std::real(Tperpm), std::imag(Tperpm));
    // printf("T_par_m(s = %.3f) = %.4e + %.4e i\n", q2, std::real(Tparm), std::imag(Tparm));
    // exit(0);

    // complex_t ALperp = A_perp(q2, -1, false);
    // complex_t ARperp = A_perp(q2, 1, false);
    // complex_t ALpar = A_par(q2, -1, false);
    // complex_t ARpar = A_par(q2, 1, false);
    // complex_t AL0 = A_0(q2, -1, false);
    // complex_t AR0 = A_0(q2, 1, false);
    // complex_t At = A_t(q2, false);
    // complex_t AS = A_S(q2, false);

    // complex_t ALperp = A_perp(q2, -1, true);
    // complex_t ARperp = A_perp(q2, 1, true);
    // complex_t ALpar = A_par(q2, -1, true);
    // complex_t ARpar = A_par(q2, 1, true);
    // complex_t AL0 = A_0(q2, -1, true);
    // complex_t AR0 = A_0(q2, 1, true);
    // complex_t At = A_t(q2, true);
    // complex_t AS = A_S(q2, true);

    // printf("A_L_perp(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ALperp), std::imag(ALperp));
	// printf("A_R_perp(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ARperp), std::imag(ARperp));
	// printf("A_L_par(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ALpar), std::imag(ALpar));
	// printf("A_R_par(s = %.3f) = %.4e + %.4e i\n", q2, std::real(ARpar), std::imag(ARpar));
	// printf("A_L_0(s = %.3f) = %.4e + %.4e i\n", q2, std::real(AL0), std::imag(AL0));
	// printf("A_R_0(s = %.3f) = %.4e + %.4e i\n", q2, std::real(AR0), std::imag(AR0));
	// printf("A_t(s = %.3f) = %.4e + %.4e i\n", q2, std::real(At), std::imag(At));
	// printf("A_S(s = %.3f) = %.4e + %.4e i\n", q2, std::real(AS), std::imag(AS));
    // exit(0);
    
    // printf("J_1c(s = %.3f) = %.4e\n", q2, J1c(q2, false));
    // printf("J_1s(s = %.3f) = %.4e\n", q2, J1s(q2, false));
    // printf("J_2c(s = %.3f) = %.4e\n", q2, J2c(q2, false));
    // printf("J_2s(s = %.3f) = %.4e\n", q2, J2s(q2, false));
    // printf("J_3(s = %.3f) = %.4e\n", q2, J3(q2, false));
    // printf("J_4(s = %.3f) = %.4e\n", q2, J4(q2, false));
    // printf("J_5(s = %.3f) = %.4e\n", q2, J5(q2, false));
    // printf("J_6c(s = %.3f) = %.4e\n", q2, J6c(q2, false));
    // printf("J_6s(s = %.3f) = %.4e\n", q2, J6s(q2, false));
    // printf("J_7(s = %.3f) = %.4e\n", q2, J7(q2, false));
    // printf("J_8(s = %.3f) = %.4e\n", q2, J8(q2, false));
    // printf("J_9(s = %.3f) = %.4e\n", q2, J9(q2, false));

    // printf("h_1c(s = %.3f) = %.4e\n", q2, h1c(q2));
    // printf("h_1s(s = %.3f) = %.4e\n", q2, h1s(q2));
    // printf("h_2c(s = %.3f) = %.4e\n", q2, h2c(q2));
    // printf("h_2s(s = %.3f) = %.4e\n", q2, h2s(q2));
    // printf("h_3(s = %.3f) = %.4e\n", q2, h3(q2));
    // printf("h_4(s = %.3f) = %.4e\n", q2, h4(q2));
    // printf("h_5(s = %.3f) = %.4e\n", q2, h5(q2));
    // printf("h_6c(s = %.3f) = %.4e\n", q2, h6c(q2));
    // printf("h_6s(s = %.3f) = %.4e\n", q2, h6s(q2));
    // printf("h_7(s = %.3f) = %.4e\n", q2, h7(q2));
    // printf("h_8(s = %.3f) = %.4e\n", q2, h8(q2));
    // printf("h_9(s = %.3f) = %.4e\n", q2, h9(q2));

    // printf("s_8(s = %.3f) = %.4e\n", q2, s8(q2));
    // printf("s_9(s = %.3f) = %.4e\n", q2, s9(q2));
    // exit(0);
}

void BsPhiDecay::fill_wilson_cache() {
    cache.C.clear();

    auto b_wilsons  = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);

    WCoef bp_cached[5] {
        WCoef::CP7, WCoef::CP9, WCoef::CP10,
        WCoef::CPQ1, WCoef::CPQ2
    };

    for (const auto& [id, val] : b_wilsons) {
        cache.C[id] = val;
    }
    for (const auto& [id, val] : bq_wilsons) {
        cache.C[id] = val;
    }
    for (auto id : bp_cached) {
        cache.C[id] = bp_wilsons.at(id);
    }
}

void BsPhiDecay::set_cfg_flags(BsPhiConfig::Lepton gen) {
    if (cfg.gen != gen) {
        cfg.gen = gen;
        load_cfg_dep_params();
    }
}

void BsPhiDecay::load_cfg_dep_params() {
    cache.m_l = (*p)(
        ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen},
        DataType::VALUE
    );

    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    cache.q2_lookup_min = std::max(cache.q2_min, 1e-4);

    auto requested_threads = cfg.n_threads;
    if (requested_threads == 0u) {
        requested_threads = std::thread::hardware_concurrency();
    }
    if (requested_threads == 0u) {
        requested_threads = 1u;
    }

    const size_t npts = BsPhiDecayCache::LOOKUP_SIZE;
    const size_t nworkers = std::min<size_t>(requested_threads, npts);

    if (nworkers <= 1u) {
        auto lam_T_perp_p = [this] (double q2, bool bar) {
            return cache.qcdf_calculator.T_perp_p(q2, bar);
        };

        fill_cache(lam_T_perp_p, cache.q2_lookup_min, cache.q2_high, cache.T_perp_p_lookup, false);
        fill_cache(lam_T_perp_p, cache.q2_lookup_min, cache.q2_high, cache.T_perp_p_bar_lookup, true);

        auto lam_T_perp_m = [this] (double q2, bool bar) {
            return cache.qcdf_calculator.T_perp_m(q2, bar);
        };

        fill_cache(lam_T_perp_m, cache.q2_lookup_min, cache.q2_high, cache.T_perp_m_lookup, false);
        fill_cache(lam_T_perp_m, cache.q2_lookup_min, cache.q2_high, cache.T_perp_m_bar_lookup, true);

        auto lam_T_par_m = [this] (double q2, bool bar) {
            return cache.qcdf_calculator.T_par_m(q2, bar);
        };

        fill_cache(lam_T_par_m, cache.q2_lookup_min, cache.q2_high, cache.T_par_m_lookup, false);
        fill_cache(lam_T_par_m, cache.q2_lookup_min, cache.q2_high, cache.T_par_m_bar_lookup, true);
    } else {
        const int B_id = 531;
        const int V_id = 333;
        const double x_min = cache.q2_lookup_min;
        const double x_max = cache.q2_high;
        const double step = (x_max - x_min) / static_cast<double>(npts - 1);

        std::vector<std::shared_ptr<BVQCDfCalculator>> qcdf_locals;
        qcdf_locals.reserve(nworkers);
        for (size_t w = 0; w < nworkers; ++w) {
            auto ff_local = std::make_shared<BVFFCalculator>(cache.ff_calculator);
            qcdf_locals.emplace_back(std::make_shared<BVQCDfCalculator>(
                B_id,
                V_id,
                cache.mu_b,
                cache.C,
                ff_local,
                cfg.ff_type,
                p,
                iobs_qcdp
            ));
        }

        std::vector<std::thread> workers;
        workers.reserve(nworkers);

        std::exception_ptr first_exception = nullptr;
        std::mutex exception_mutex;

        auto worker = [&] (size_t worker_id, size_t begin, size_t end) {
            try {
                BVQCDfCalculator& qcdf = *qcdf_locals[worker_id];
                for (size_t i = begin; i < end; ++i) {
                    const double q2 = x_min + step * static_cast<double>(i);

                    cache.T_perp_p_lookup[i]     = qcdf.T_perp_p(q2, false);
                    cache.T_perp_p_bar_lookup[i] = qcdf.T_perp_p(q2, true);

                    cache.T_perp_m_lookup[i]     = qcdf.T_perp_m(q2, false);
                    cache.T_perp_m_bar_lookup[i] = qcdf.T_perp_m(q2, true);

                    cache.T_par_m_lookup[i]      = qcdf.T_par_m(q2, false);
                    cache.T_par_m_bar_lookup[i]  = qcdf.T_par_m(q2, true);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(exception_mutex);
                if (!first_exception) {
                    first_exception = std::current_exception();
                }
            }
        };

        const size_t chunk = (npts + nworkers - 1) / nworkers;
        for (size_t w = 0; w < nworkers; ++w) {
            const size_t begin = w * chunk;
            const size_t end = std::min(npts, begin + chunk);
            if (begin >= end) {
                break;
            }
            workers.emplace_back(worker, w, begin, end);
        }

        for (auto& th : workers) {
            th.join();
        }

        if (first_exception) {
            std::rethrow_exception(first_exception);
        }
    }

    compute_binned_J_i();
}

complex_t BsPhiDecay::T_perp_p_cached(double q2, bool bar) {
    const double x = std::max(q2, cache.q2_lookup_min);
    return lerp(
        x,
        bar ? cache.T_perp_p_bar_lookup : cache.T_perp_p_lookup,
        cache.q2_lookup_min,
        cache.q2_high
    );
}

complex_t BsPhiDecay::T_perp_m_cached(double q2, bool bar) {
    const double x = std::max(q2, cache.q2_lookup_min);
    return lerp(
        x,
        bar ? cache.T_perp_m_bar_lookup : cache.T_perp_m_lookup,
        cache.q2_lookup_min,
        cache.q2_high
    );
}

complex_t BsPhiDecay::T_par_m_cached(double q2, bool bar) {
    const double x = std::max(q2, cache.q2_lookup_min);
    return lerp(
        x,
        bar ? cache.T_par_m_bar_lookup : cache.T_par_m_lookup,
        cache.q2_lookup_min,
        cache.q2_high
    );
}

double BsPhiDecay::beta_l(double q2) {
    const double x = 1.0 - 4.0 * cache.m_l * cache.m_l / q2;
    return std::sqrt(std::max(0.0, x));
}

double BsPhiDecay::lambda(double q2) {
    const double mB2 = cache.m_Bs * cache.m_Bs;
    const double mphi2 = cache.m_phi * cache.m_phi;

    const double lam =
        mB2 * mB2
        + mphi2 * mphi2
        + q2 * q2
        - 2.0 * (mB2 * mphi2 + (mB2 + mphi2) * q2);

    return std::max(0.0, lam);
}

complex_t BsPhiDecay::N(double q2, bool bar) {
    complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0; 
    return N0 * std::sqrt(q2 * beta_l(q2) * std::sqrt(lambda(q2)));
}

complex_t BsPhiDecay::delta_A_perp(double q2, double sign, bool bar) {
    size_t id = size_t (0.5 * (1 + sign));
    complex_t guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    double m_b = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    return 2 * RT2 * m_b * N(q2, bar) * std::sqrt(lambda(q2)) / q2 * T_perp_p_cached(q2, bar) * guesstimate_err;
}

complex_t BsPhiDecay::delta_A_par(double q2, double sign, bool bar) {
    size_t id = 2 + size_t (0.5 * (1 + sign));
    complex_t guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    double m_b = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    return -4 * RT2 * m_b * N(q2, bar) * (cache.m_Bs * cache.m_Bs - cache.m_phi * cache.m_phi) * cache.ff_calculator.E(q2) / (q2 * cache.m_Bs) * T_perp_m_cached(q2, bar) * guesstimate_err;
}

complex_t BsPhiDecay::A_perp_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_Bs + cache.m_phi);
        F_T = T_perp_p_cached(q2, bar);
        size_t id = size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BV_FF::T1, q2);
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        w += cache.qcdf_calculator.Y(q2);
        F = w * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_Bs + cache.m_phi);
        delta_A = delta_A_perp(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("----------- A_%c_perp_low(s = %.3f) -----------\n", sign == -1 ? 'L' : 'R', q2);
    // printf("prefactor = %.4e + %.4e i\n", real(N(q2, bar) * std::sqrt(2 * lambda(q2))), imag(N(q2, bar) * std::sqrt(2 * lambda(q2))));
    // printf("F = %.4e + %.4e i\n", real(F), imag(F));
    // printf("m_b_local = %.4e\n", m_b_local);
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));
    // printf("delta_A = %.4e + %.4e i\n", real(delta_A), imag(delta_A));

    return (N(q2, bar) * std::sqrt(2 * lambda(q2)) * (F + 2. * m_b_local * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BsPhiDecay::A_par_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * 2. * cache.ff_calculator.E(q2) / (cache.m_Bs + cache.m_phi) * cache.ff_calculator.get(BV_FF::XI_PERP, q2) / (cache.m_Bs - cache.m_phi);
        F_T = 2. * cache.ff_calculator.E(q2) * T_perp_m_cached(q2, bar) / cache.m_Bs;
        size_t id = 2 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * cache.ff_calculator.get(BV_FF::T2, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        w += cache.qcdf_calculator.Y(q2);
        F = w * cache.ff_calculator.get(BV_FF::A1, q2) / (cache.m_Bs - cache.m_phi);
        delta_A = delta_A_par(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("----------- A_%c_par_low(s = %.3f) -----------\n", sign == -1 ? 'L' : 'R', q2);
    // printf("prefactor = %.4e + %.4e i\n", real(-N(q2, bar) * std::sqrt(2.) * (cache.m_Bs * cache.m_Bs - cache.m_phi * cache.m_phi)), imag(-N(q2, bar) * std::sqrt(2.) * (cache.m_Bs * cache.m_Bs - cache.m_phi * cache.m_phi)));
    // printf("F = %.4e + %.4e i\n", real(F), imag(F));
    // printf("m_b_local = %.4e\n", m_b_local);
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));
    // printf("delta_A = %.4e + %.4e i\n", real(delta_A), imag(delta_A));

    return (-N(q2, bar) * std::sqrt(2.) * (cache.m_Bs * cache.m_Bs - cache.m_phi * cache.m_phi) * (F + 2. * m_b_local * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BsPhiDecay::A_0_low(double q2, double sign, bool bar) {
    double mB2 = cache.m_Bs * cache.m_Bs;
    double mK2 = cache.m_phi * cache.m_phi;
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * ((mB2 - mK2 - q2) * 2. * cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PERP, q2) - lambda(q2) * cache.m_Bs / (mB2 - mK2) * (cache.ff_calculator.get(BV_FF::XI_PERP, q2) - cache.ff_calculator.get(BV_FF::XI_PAR, q2)));
        F_T = 2. * cache.ff_calculator.E(q2) * (mB2 + 3. * mK2 - q2) * T_perp_m_cached(q2, bar) / cache.m_Bs - lambda(q2) * (T_perp_m_cached(q2, bar) + T_par_m_cached(q2, bar)) / (mB2 - mK2);
        size_t id = 4 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * 8. * cache.m_Bs * mK2 / (cache.m_Bs + cache.m_phi) * cache.ff_calculator.get(BV_FF::T23, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        w += cache.qcdf_calculator.Y(q2);
        F = w * 16. * cache.m_Bs * mK2 * cache.ff_calculator.get(BV_FF::A12, q2);
        delta_A = delta_A_0(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("----------- A_%c_0_low(s = %.3f) -----------\n", sign == -1 ? 'L' : 'R', q2);
    // printf("prefactor = %.4e + %.4e i\n", real(-N(q2, bar) / (2. * cache.m_phi * std::sqrt(q2))), imag(-N(q2, bar) / (2. * cache.m_phi * std::sqrt(q2))));
    // printf("F = %.4e + %.4e i\n", real(F), imag(F));
    // printf("m_b_local = %.4e\n", m_b_local);
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));
    // printf("delta_A = %.4e + %.4e i\n", real(delta_A), imag(delta_A));

    return (-N(q2, bar) / (2. * cache.m_phi * std::sqrt(q2)) * (F + 2. * m_b_local * F_T) + delta_A) * had_err_factor;
}

complex_t BsPhiDecay::A_t_low(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }

    complex_t F;
    if (cfg.ff_type == B_FF_Type::SOFT) {
        F = cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PAR, q2) / (cache.m_phi * cache.qcdf_calculator.Delta_par(q2));
    } else {
        F = cache.ff_calculator.get(BV_FF::A0, q2);
    }

    return N(q2, bar) * std::sqrt(lambda(q2) / q2) * ((2. * C10 + q2 / (cache.m_l * (cache.m_b_mu_b + cache.m_s)) * CQ2)) * F;
}

complex_t BsPhiDecay::A_S_low(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);

    complex_t F;
    if (cfg.ff_type == B_FF_Type::SOFT) {
        F = cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PAR, q2) / (cache.m_phi * cache.qcdf_calculator.Delta_par(q2));
    } else {
        F = cache.ff_calculator.get(BV_FF::A0, q2);
    }

    return -2. * N(q2, bar) * std::sqrt(lambda(q2)) * CQ1 / (cache.m_b_mu_b + cache.m_s) * F;
}

complex_t BsPhiDecay::delta_A_0(double q2, double sign, bool bar) {
    size_t id = 4 + size_t (0.5 * (1 + sign));
    complex_t guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;

    double mB2 = cache.m_Bs * cache.m_Bs;
    double mB3 = cache.m_Bs * mB2;
    double mK2 = cache.m_phi * cache.m_phi;
    double f = lambda(q2) / ((mB2 - mK2) * mB2);
    double m_b = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    complex_t delta_A = -N(q2, bar) * m_b * mB2 / (std::sqrt(q2) * cache.m_phi) * (
        (2 * (mB2 + 3 * mK2 - q2) * cache.ff_calculator.E(q2) / mB3 - f) * T_perp_m_cached(q2, bar)
      - f * T_par_m_cached(q2, bar)
    );
    return delta_A * guesstimate_err;
}

complex_t BsPhiDecay::C7_eff(double q2, bool bar) {
    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    return (!bar ? std::conj(cache.C[WCoef::C7]) : cache.C[WCoef::C7]) + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*BV::f_87(s_hat, cache.L_b));
}

complex_t BsPhiDecay::C9_eff(double q2, bool bar) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+l_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t B = BV::B_Seidel(s_hat, cache.L_b);
    complex_t C = BV::C_Seidel(q2, cache.mu_b);

    return (!bar ? std::conj(cache.C[WCoef::C9]) : cache.C[WCoef::C9])
         + BV::h(q2, 0., cache.mu_b) * C_h0
         + BV::h(q2, cache.m_b_PS, cache.mu_b) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * BV::f_89(s_hat))
         + std::pow(cache.m_c_m_c, 2) / q2 * C_mc;
}

complex_t BsPhiDecay::A_perp_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) + (!bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) + (!bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] + cache.C[WCoef::CP10];
    if (!bar) C10 = std::conj(C10);

    // printf("C7_eff = %.4e + %.4e i\n", real(C7), imag(C7));
    // printf("C9_eff = %.4e + %.4e i\n", real(C9), imag(C9));
    // printf("C10_eff = %.4e + %.4e i\n", real(C10), imag(C10));
    // printf("f_perp = %.4e\n", cache.ff_calculator.get(BV_FF::F_PERP, q2));
    // printf("N = %.4e + %.4e i\n", real(N(q2, bar)), imag(N(q2, bar)));
    // printf("pref T = %.4e\n", 2. * cache.kappa * cache.m_b_m_b * cache.m_Bs / q2);

    return N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_m_b * cache.m_Bs / q2 * C7) * cache.ff_calculator.get(BV_FF::F_PERP, q2) * (1. + cache.A_had_err_high[size_t (0.5 * (1 + sign))]);
}

complex_t BsPhiDecay::A_par_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (!bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (!bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (!bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_m_b * cache.m_Bs / q2 * C7) * cache.ff_calculator.get(BV_FF::F_PAR, q2) * (1. + cache.A_had_err_high[2 + size_t (0.5 * (1 + sign))]);
}

complex_t BsPhiDecay::A_0_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (!bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (!bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (!bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_m_b * cache.m_Bs / q2 * C7) * cache.ff_calculator.get(BV_FF::F_0, q2)  * (1. + cache.A_had_err_high[4 + size_t (0.5 * (1 + sign))]);
}

complex_t BsPhiDecay::A_t_high(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (!bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }
    return N(q2, bar) * sqrt(lambda(q2) / q2) * (2. * C10 + q2 / cache.m_l * CQ2 / (cache.m_b_m_b + cache.m_s)) * cache.ff_calculator.get(BV_FF::A0, q2) * (1. + cache.A_had_err_high[6]);
}

complex_t BsPhiDecay::A_S_high(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (!bar) CQ1 = std::conj(CQ1);
    return -2. * N(q2, bar) * sqrt(lambda(q2)) * CQ1 / (cache.m_b_m_b + cache.m_s) * cache.ff_calculator.get(BV_FF::A0, q2) * (1. + cache.A_had_err_high[7]);
}

complex_t BsPhiDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BsPhiDecay::A_perp(double q2, double sign, bool bar) {
    return interpolate(q2, A_perp_low(q2, sign, !bar), A_perp_high(q2, sign, !bar));
}

complex_t BsPhiDecay::A_par(double q2, double sign, bool bar) {
    return interpolate(q2, A_par_low(q2, sign, !bar), A_par_high(q2, sign, !bar));
}

complex_t BsPhiDecay::A_0(double q2, double sign, bool bar) {
    return interpolate(q2, A_0_low(q2, sign, !bar), A_0_high(q2, sign, !bar));
}

complex_t BsPhiDecay::A_t(double q2, bool bar) {
    return interpolate(q2, A_t_low(q2, !bar), A_t_high(q2, !bar));
}

complex_t BsPhiDecay::A_S(double q2, bool bar) {
    return interpolate(q2, A_S_low(q2, !bar), A_S_high(q2, !bar));
}

double BsPhiDecay::J1s(double q2, bool bar) {
    return (2. + std::pow(beta_l(q2), 2)) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_par(q2, -1, bar)), 2)
      + std::pow(std::abs(A_par(q2, 1, bar)), 2)
    ) + std::pow(2. * cache.m_l, 2) / q2 * std::real(
        A_perp(q2, -1, bar) * std::conj(A_perp(q2, 1, bar))
      + A_par(q2, -1, bar) * std::conj(A_par(q2, 1, bar))
    );
}

double BsPhiDecay::J1c(double q2, bool bar) {
    return std::pow(std::abs(A_0(q2, -1, bar)), 2) + std::pow(std::abs(A_0(q2, 1, bar)), 2) 
         + std::pow(2 * cache.m_l, 2) / q2 * (
              std::pow(std::abs(A_t(q2, bar)), 2)
            + 2. * std::real(A_0(q2, -1, bar) * std::conj(A_0(q2, 1, bar)))
           )
         + std::pow(beta_l(q2) * std::abs(A_S(q2, bar)), 2);
}

double BsPhiDecay::J2s(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_par(q2, -1, bar)), 2)
      + std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BsPhiDecay::J2c(double q2, bool bar) {
    return -std::pow(beta_l(q2), 2) * (
        std::pow(std::abs(A_0(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_0(q2, 1, bar)), 2)
    );
}

double BsPhiDecay::J3(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 2. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      - std::pow(std::abs(A_par(q2, -1, bar)), 2)
      - std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BsPhiDecay::J4(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::real(A_0(q2, 1, bar) * std::conj(A_par(q2, 1, bar)))
    );
}

double BsPhiDecay::J5(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
      - cache.m_l / std::sqrt(q2) * std::real((A_par(q2, -1, bar) + A_par(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BsPhiDecay::J6s(double q2, bool bar) {
    return 2. * beta_l(q2) * (
        std::real(A_par(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_par(q2, 1, bar) * std::conj(A_perp(q2, 1, bar))) 
    );
}

double BsPhiDecay::J6c(double q2, bool bar) {
    return 4. * beta_l(q2) * cache.m_l / std::sqrt(q2) * (std::real((A_0(q2, -1, bar) + A_0(q2, 1, bar)) * std::conj(A_S(q2, bar))));
}

double BsPhiDecay::J7(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      - std::imag(A_0(q2, 1, bar) * std::conj(A_par(q2, 1, bar)))
      + cache.m_l / std::sqrt(q2) * std::imag((A_perp(q2, -1, bar) + A_perp(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BsPhiDecay::J8(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      + std::imag(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
    );
}

double BsPhiDecay::J9(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) * (
        std::imag(A_perp(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::imag(A_perp(q2, 1, bar) * std::conj(A_par(q2, 1, bar))) 
    );
}

double BsPhiDecay::h1s(double q2) {
    complex_t ALperp = A_perp(q2, -1, false);
    complex_t ARperp = A_perp(q2, 1, false);
    complex_t ALpar = A_par(q2, -1, false);
    complex_t ARpar = A_par(q2, 1, false);
    complex_t ALperp_tilde = -A_perp(q2, -1, true);
    complex_t ARperp_tilde = -A_perp(q2, 1, true);
    complex_t ALpar_tilde = A_par(q2, -1, true);
    complex_t ARpar_tilde = A_par(q2, 1, true);
    complex_t sh1s_A = cache.up * (ALperp_tilde * std::conj(ALperp) + ALpar_tilde * std::conj(ALpar) + ARperp_tilde * std::conj(ARperp) + ARpar_tilde * std::conj(ARpar));
	complex_t sh1s_B = cache.up * (ALperp_tilde * std::conj(ARperp) + ALpar_tilde * std::conj(ARpar));
	complex_t sh1s_C = cache.um * (ALperp * std::conj(ARperp_tilde) + ALpar * std::conj(ARpar_tilde));
    double x = std::pow(cache.m_l, 2) / q2;

    return (2. + std::pow(beta_l(q2), 2)) / 2. * std::real(sh1s_A) + 4. * x * std::real(sh1s_B + sh1s_C);
}

double BsPhiDecay::h1c(double q2) {
    complex_t AL0 = A_0(q2, -1, false);
    complex_t AR0 = A_0(q2, 1, false);
    complex_t AL0_tilde = A_0(q2, -1, true);
    complex_t AR0_tilde = A_0(q2, 1, true);
    complex_t sh1c_A = cache.up * (AL0_tilde * std::conj(AL0) + AR0_tilde * std::conj(AR0));
	complex_t sh1c_B = cache.up * A_t(q2, true) * std::conj(A_t(q2, false));
	complex_t sh1c_C = cache.up * AL0_tilde * std::conj(AR0);
	complex_t sh1c_D = cache.um * AL0 * std::conj(AR0_tilde);
	complex_t sh1c_E = cache.up * -A_S(q2, true) * std::conj(A_S(q2, false));
    double x = std::pow(cache.m_l, 2) / q2;

    return 2. * std::real(sh1c_A) + 8. * x * (std::real(sh1c_B) + std::real(sh1c_C + sh1c_D)) + 2. * std::pow(beta_l(q2), 2) * std::real(sh1c_E);
}

double BsPhiDecay::h2s(double q2) {
    complex_t sh1s_A = cache.up * (-A_perp(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, -1, true) * std::conj(A_par(q2, -1, false)) - A_perp(q2, 1, true) * std::conj(A_perp(q2, 1, false)) + A_par(q2, 1, true) * std::conj(A_par(q2, 1, false)));
    return std::pow(beta_l(q2), 2) / 2. * std::real(sh1s_A);
}

double BsPhiDecay::h2c(double q2) {
    complex_t sh1c_A = cache.up * (A_0(q2, -1, true) * std::conj(A_0(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_0(q2, 1, false)));
    return -2. * std::pow(beta_l(q2), 2) * std::real(sh1c_A);
}

double BsPhiDecay::h3(double q2) {
    complex_t sh3_A = -cache.up * (A_perp(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, -1, true) * std::conj(A_par(q2, -1, false)) + A_perp(q2, 1, true) * std::conj(A_perp(q2, 1, false)) + A_par(q2, 1, true) * std::conj(A_par(q2, 1, false)));
    return std::pow(beta_l(q2), 2) * std::real(sh3_A);
}

double BsPhiDecay::h4(double q2) {
    complex_t sh4_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_par(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_par(q2, 1, false)));
	complex_t sh4_B  = cache.um * (A_0(q2, -1, false) * std::conj(A_par(q2, -1, true)) + A_0(q2, 1, false) * std::conj(A_par(q2, 1, true)));
    return INV_RT2 * std::pow(beta_l(q2), 2) * std::real(sh4_A + sh4_B);
}

double BsPhiDecay::h5(double q2) {
    complex_t sh5_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_perp(q2, -1, false)) - A_0(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh5_B  = cache.um * (A_0(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) - A_0(q2, 1, false) * std::conj(A_par(q2, 1, true)));
	complex_t sh5_C  = cache.up * ((A_par(q2, -1, true) + A_par(q2, 1, true)) * std::conj(A_S(q2, false)));
	complex_t sh5_D  = cache.um * ((A_par(q2, -1, false) + A_par(q2, 1, false)) * std::conj(-A_S(q2, true)));
    return RT2 * beta_l(q2) * (std::real(sh5_A + sh5_B)  - cache.m_l / std::sqrt(q2) * std::real(sh5_C + sh5_D));
}

double BsPhiDecay::h6s(double q2) {
    complex_t sh6s_A = cache.up * (A_par(q2, -1, true) * std::conj(A_perp(q2, -1, false)) - A_par(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh6s_B = cache.um * (A_par(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) - A_par(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return 2. * beta_l(q2) * std::real(sh6s_A + sh6s_B);
}

double BsPhiDecay::h6c(double q2) {
    complex_t sh6c_A = cache.up * ((A_0(q2, -1, true) + A_0(q2, 1, true)) * std::conj(A_S(q2, false)));
	complex_t sh6c_B = cache.um * ((A_0(q2, -1, false) + A_0(q2, 1, false)) * std::conj(-A_S(q2, true)));
    return 4. * beta_l(q2) * cache.m_l / std::sqrt(q2) * std::real(sh6c_A + sh6c_B);
}

double BsPhiDecay::h7(double q2) {
    complex_t sh7_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_par(q2, -1, false)) - A_0(q2, 1, true) * std::conj(A_par(q2, 1, false)));
	complex_t sh7_B  = cache.um * (A_0(q2, -1, false) * std::conj(A_par(q2, -1, true)) - A_0(q2, 1, false) * std::conj(A_par(q2, 1, true)));
	complex_t sh7_C  = cache.up * (-(A_perp(q2, -1, true) + A_perp(q2, 1, true)) * std::conj(A_S(q2, false)));
	complex_t sh7_D  = cache.um * ((A_perp(q2, -1, false) + A_perp(q2, 1, false)) * std::conj(-A_S(q2, true)));
    return RT2 * beta_l(q2) * (std::imag(sh7_A + sh7_B) + cache.m_l / std::sqrt(q2) * std::imag(sh7_C + sh7_D));
}

double BsPhiDecay::h8(double q2) {
    complex_t sh8_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh8_B  = cache.um * (A_0(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_0(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return INV_RT2 * std::pow(beta_l(q2), 2) * std::imag(sh8_A + sh8_B);
}

double BsPhiDecay::h9(double q2) {
    complex_t sh9_A  = cache.up * (A_par(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh9_B  = cache.um * (A_par(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_par(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return -std::pow(beta_l(q2), 2) * std::imag(sh9_A + sh9_B);
}

double BsPhiDecay::s8(double q2) {
    complex_t sh8_A  = cache.up * (A_0(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_0(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh8_B  = cache.um * (A_0(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_0(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return -1. / sqrt(q2) * std::pow(beta_l(q2), 2) * std::real(sh8_A - sh8_B);
}

double BsPhiDecay::s9(double q2) {
    complex_t sh9_A  = cache.up * (A_par(q2, -1, true) * std::conj(A_perp(q2, -1, false)) + A_par(q2, 1, true) * std::conj(A_perp(q2, 1, false)));
	complex_t sh9_B  = cache.um * (A_par(q2, -1, false) * std::conj(-A_perp(q2, -1, true)) + A_par(q2, 1, false) * std::conj(-A_perp(q2, 1, true)));
    return std::pow(beta_l(q2), 2) * std::real(sh9_A - sh9_B);
}

void BsPhiDecay::compute_binned_J_i() {
    // One pass per bin:
    // - evaluate all transversity amplitudes once for the direct mode (bar=false)
    //   and once for the CP-conjugate mode (bar=true)
    // - derive J_i, h_i and s_i from those amplitudes
    // - accumulate every needed binned integral together with a fixed GL24 quadrature
    //
    // This removes the previous pattern of 15 independent adaptive integrals per bin,
    // each of which was recomputing the same amplitudes many times.
    static constexpr std::array<double, 24> GL24_X {{
        -0.99518721999702131, -0.97472855597130947, -0.93827455200273280, -0.88641552700440107,
        -0.82000198597390295, -0.74012419157855436, -0.64809365193697555, -0.54542147138883956,
        -0.43379350762604513, -0.31504267969616340, -0.19111886747361631, -0.06405689286260563,
         0.06405689286260563,  0.19111886747361631,  0.31504267969616340,  0.43379350762604513,
         0.54542147138883956,  0.64809365193697555,  0.74012419157855436,  0.82000198597390295,
         0.88641552700440107,  0.93827455200273280,  0.97472855597130947,  0.99518721999702131
    }};
    static constexpr std::array<double, 24> GL24_W {{
        0.01234122979998869, 0.02853138862893356, 0.04427743881741941, 0.05929858491543636,
        0.07334648141108016, 0.08619016153195321, 0.09761865210411393, 0.10744427011596556,
        0.11550566805372552, 0.12167047292780329, 0.12583745634682825, 0.12793819534675202,
        0.12793819534675202, 0.12583745634682825, 0.12167047292780329, 0.11550566805372552,
        0.10744427011596556, 0.09761865210411393, 0.08619016153195321, 0.07334648141108016,
        0.05929858491543636, 0.04427743881741941, 0.02853138862893356, 0.01234122979998869
    }};

    struct AmpSet {
        complex_t A_perp_L;
        complex_t A_perp_R;
        complex_t A_par_L;
        complex_t A_par_R;
        complex_t A_0_L;
        complex_t A_0_R;
        complex_t A_t;
        complex_t A_S;
        double beta;
        double beta2;
        double sqrt_q2;
    };

    struct JSet {
        double j1s;
        double j1c;
        double j2s;
        double j2c;
        double j3;
        double j4;
        double j5;
        double j6s;
        double j6c;
        double j7;
        double j8;
        double j9;
    };

    auto zero_if_close = [] (double x, double tol) {
        return std::abs(x) < tol ? 0.0 : x;
    };

    auto clear_and_reserve = [&] (std::array<std::vector<double>, 15>& dest, size_t nbins) {
        for (auto& v : dest) {
            v.clear();
            v.reserve(nbins);
        }
    };

    auto eval_amplitudes = [&] (double q2, bool bar) -> AmpSet {
        const double beta = beta_l(q2);
        return {
            A_perp(q2, -1, bar),
            A_perp(q2,  1, bar),
            A_par(q2,  -1, bar),
            A_par(q2,   1, bar),
            A_0(q2,    -1, bar),
            A_0(q2,     1, bar),
            A_t(q2, bar),
            A_S(q2, bar),
            beta,
            beta * beta,
            std::sqrt(q2)
        };
    };

    auto eval_J = [&] (const AmpSet& a, double q2) -> JSet {
        const double sqrt2 = std::sqrt(2.0);
        const double ml = cache.m_l;
        const double ml_over_sqrt_q2 = ml / a.sqrt_q2;
        const double four_ml2_over_q2 = 4.0 * ml * ml / q2;

        const double norm_Aperp_L = std::norm(a.A_perp_L);
        const double norm_Aperp_R = std::norm(a.A_perp_R);
        const double norm_Apar_L  = std::norm(a.A_par_L);
        const double norm_Apar_R  = std::norm(a.A_par_R);
        const double norm_A0_L    = std::norm(a.A_0_L);
        const double norm_A0_R    = std::norm(a.A_0_R);
        const double norm_At      = std::norm(a.A_t);
        const double norm_AS      = std::norm(a.A_S);

        const complex_t Apar_sum  = a.A_par_L + a.A_par_R;
        const complex_t A0_sum    = a.A_0_L + a.A_0_R;
        const complex_t Aperp_sum = a.A_perp_L + a.A_perp_R;

        JSet j {};
        j.j1s =
            (2.0 + a.beta2) / 4.0 * (norm_Aperp_L + norm_Aperp_R + norm_Apar_L + norm_Apar_R)
            + four_ml2_over_q2 * std::real(a.A_perp_L * std::conj(a.A_perp_R)
                                         + a.A_par_L  * std::conj(a.A_par_R));

        j.j1c =
            norm_A0_L + norm_A0_R
            + four_ml2_over_q2 * (norm_At + 2.0 * std::real(a.A_0_L * std::conj(a.A_0_R)))
            + a.beta2 * norm_AS;

        j.j2s = a.beta2 / 4.0 * (norm_Aperp_L + norm_Aperp_R + norm_Apar_L + norm_Apar_R);
        j.j2c = -a.beta2 * (norm_A0_L + norm_A0_R);

        j.j3 = a.beta2 / 2.0 * (norm_Aperp_L + norm_Aperp_R - norm_Apar_L - norm_Apar_R);

        j.j4 = a.beta2 / sqrt2 * (
            std::real(a.A_0_L * std::conj(a.A_par_L))
          + std::real(a.A_0_R * std::conj(a.A_par_R))
        );

        j.j5 = a.beta * sqrt2 * (
            std::real(a.A_0_L * std::conj(a.A_perp_L))
          - std::real(a.A_0_R * std::conj(a.A_perp_R))
          - ml_over_sqrt_q2 * std::real(Apar_sum * std::conj(a.A_S))
        );

        j.j6s = 2.0 * a.beta * (
            std::real(a.A_par_L * std::conj(a.A_perp_L))
          - std::real(a.A_par_R * std::conj(a.A_perp_R))
        );

        j.j6c = 4.0 * a.beta * ml_over_sqrt_q2 * std::real(A0_sum * std::conj(a.A_S));

        j.j7 = a.beta * sqrt2 * (
            std::imag(a.A_0_L * std::conj(a.A_par_L))
          - std::imag(a.A_0_R * std::conj(a.A_par_R))
          + ml_over_sqrt_q2 * std::imag(Aperp_sum * std::conj(a.A_S))
        );

        j.j8 = a.beta2 / sqrt2 * (
            std::imag(a.A_0_L * std::conj(a.A_perp_L))
          + std::imag(a.A_0_R * std::conj(a.A_perp_R))
        );

        j.j9 = a.beta2 * (
            std::imag(a.A_perp_L * std::conj(a.A_par_L))
          + std::imag(a.A_perp_R * std::conj(a.A_par_R))
        );

        return j;
    };

    auto eval_integrands = [&] (double q2) -> std::array<double, 15> {
        const AmpSet a = eval_amplitudes(q2, false);
        const AmpSet b = eval_amplitudes(q2, true);
        const JSet jf = eval_J(a, q2);
        const JSet jt = eval_J(b, q2);

        const double ml = cache.m_l;
        const double x = ml * ml / q2;
        const double beta = a.beta;
        const double beta2 = a.beta2;
        const double sqrt_q2 = a.sqrt_q2;

        const complex_t ALperp_tilde = -b.A_perp_L;
        const complex_t ARperp_tilde = -b.A_perp_R;
        const complex_t ALpar_tilde  =  b.A_par_L;
        const complex_t ARpar_tilde  =  b.A_par_R;
        const complex_t AL0_tilde    =  b.A_0_L;
        const complex_t AR0_tilde    =  b.A_0_R;
        const complex_t At_tilde     =  b.A_t;
        const complex_t AS_tilde     = -b.A_S;

        const complex_t sh1s_A = cache.up * (ALperp_tilde * std::conj(a.A_perp_L) + ALpar_tilde * std::conj(a.A_par_L)
                                           + ARperp_tilde * std::conj(a.A_perp_R) + ARpar_tilde * std::conj(a.A_par_R));
        const complex_t sh1s_B = cache.up * (ALperp_tilde * std::conj(a.A_perp_R) + ALpar_tilde * std::conj(a.A_par_R));
        const complex_t sh1s_C = cache.um * (a.A_perp_L * std::conj(ARperp_tilde) + a.A_par_L * std::conj(ARpar_tilde));
        const double h1s = (2.0 + beta2) / 2.0 * std::real(sh1s_A) + 4.0 * x * std::real(sh1s_B + sh1s_C);

        const complex_t sh1c_A = cache.up * (AL0_tilde * std::conj(a.A_0_L) + AR0_tilde * std::conj(a.A_0_R));
        const complex_t sh1c_B = cache.up * At_tilde * std::conj(a.A_t);
        const complex_t sh1c_C = cache.up * AL0_tilde * std::conj(a.A_0_R);
        const complex_t sh1c_D = cache.um * a.A_0_L * std::conj(AR0_tilde);
        const complex_t sh1c_E = cache.up * AS_tilde * std::conj(a.A_S);
        const double h1c = 2.0 * std::real(sh1c_A)
                         + 8.0 * x * (std::real(sh1c_B) + std::real(sh1c_C + sh1c_D))
                         + 2.0 * beta2 * std::real(sh1c_E);

        const complex_t sh2s_A = cache.up * (ALperp_tilde * std::conj(a.A_perp_L) + ALpar_tilde * std::conj(a.A_par_L)
                                           + ARperp_tilde * std::conj(a.A_perp_R) + ARpar_tilde * std::conj(a.A_par_R));
        const double h2s = beta2 / 2.0 * std::real(sh2s_A);

        const complex_t sh2c_A = cache.up * (AL0_tilde * std::conj(a.A_0_L) + AR0_tilde * std::conj(a.A_0_R));
        const double h2c = -2.0 * beta2 * std::real(sh2c_A);

        const complex_t sh3_A = -cache.up * (b.A_perp_L * std::conj(a.A_perp_L) + b.A_par_L * std::conj(a.A_par_L)
                                           + b.A_perp_R * std::conj(a.A_perp_R) + b.A_par_R * std::conj(a.A_par_R));
        const double h3 = beta2 * std::real(sh3_A);

        const complex_t sh4_A = cache.up * (b.A_0_L * std::conj(a.A_par_L) + b.A_0_R * std::conj(a.A_par_R));
        const complex_t sh4_B = cache.um * (a.A_0_L * std::conj(b.A_par_L) + a.A_0_R * std::conj(b.A_par_R));
        const double h4 = INV_RT2 * beta2 * std::real(sh4_A + sh4_B);

        const complex_t sh5_A = cache.up * (b.A_0_L * std::conj(a.A_perp_L) - b.A_0_R * std::conj(a.A_perp_R));
        const complex_t sh5_B = cache.um * (a.A_0_L * std::conj(-b.A_perp_L) - a.A_0_R * std::conj(b.A_par_R));
        const complex_t sh5_C = cache.up * ((b.A_par_L + b.A_par_R) * std::conj(a.A_S));
        const complex_t sh5_D = cache.um * ((a.A_par_L + a.A_par_R) * std::conj(-b.A_S));
        const double h5 = RT2 * beta * (std::real(sh5_A + sh5_B) - ml / sqrt_q2 * std::real(sh5_C + sh5_D));

        const complex_t sh6s_A = cache.up * (b.A_par_L * std::conj(a.A_perp_L) - b.A_par_R * std::conj(a.A_perp_R));
        const complex_t sh6s_B = cache.um * (a.A_par_L * std::conj(-b.A_perp_L) - a.A_par_R * std::conj(-b.A_perp_R));
        const double h6s = 2.0 * beta * std::real(sh6s_A + sh6s_B);

        const complex_t sh6c_A = cache.up * ((b.A_0_L + b.A_0_R) * std::conj(a.A_S));
        const complex_t sh6c_B = cache.um * ((a.A_0_L + a.A_0_R) * std::conj(-b.A_S));
        const double h6c = 4.0 * beta * ml / sqrt_q2 * std::real(sh6c_A + sh6c_B);

        const complex_t sh7_A = cache.up * (b.A_0_L * std::conj(a.A_par_L) - b.A_0_R * std::conj(a.A_par_R));
        const complex_t sh7_B = cache.um * (a.A_0_L * std::conj(b.A_par_L) - a.A_0_R * std::conj(b.A_par_R));
        const complex_t sh7_C = cache.up * (-(b.A_perp_L + b.A_perp_R) * std::conj(a.A_S));
        const complex_t sh7_D = cache.um * ((a.A_perp_L + a.A_perp_R) * std::conj(-b.A_S));
        const double h7 = RT2 * beta * (std::imag(sh7_A + sh7_B) + ml / sqrt_q2 * std::imag(sh7_C + sh7_D));

        const complex_t sh8_A = cache.up * (b.A_0_L * std::conj(a.A_perp_L) + b.A_0_R * std::conj(a.A_perp_R));
        const complex_t sh8_B = cache.um * (a.A_0_L * std::conj(-b.A_perp_L) + a.A_0_R * std::conj(-b.A_perp_R));
        const double h8 = INV_RT2 * beta2 * std::imag(sh8_A + sh8_B);

        const complex_t sh9_A = cache.up * (b.A_par_L * std::conj(a.A_perp_L) + b.A_par_R * std::conj(a.A_perp_R));
        const complex_t sh9_B = cache.um * (a.A_par_L * std::conj(-b.A_perp_L) + a.A_par_R * std::conj(-b.A_perp_R));
        const double h9 = -beta2 * std::imag(sh9_A + sh9_B);

        const double s8v = -beta2 / sqrt_q2 * std::real(sh8_A - sh8_B);
        const double s9v =  beta2 * std::real(sh9_A - sh9_B);

        const double val6  = zero_if_close(jf.j5  - jt.j5  - cache.ys * h5,  1e-30);
        const double val7  = zero_if_close(jf.j6s - jt.j6s - cache.ys * h6s, 1e-30);
        const double val9  = zero_if_close(jf.j6c - jt.j6c - cache.ys * h6c, 1e-30);
        const double val11 = zero_if_close(jf.j8  - jt.j8  - cache.ys * h8,  1e-30);
        const double val12 = zero_if_close(jf.j9  - jt.j9  - cache.ys * h9,  1e-30);

        return {{
            2.0 * (jf.j1s - cache.ys * h1s / 2.0) + jf.j1c - cache.ys * h1c / 2.0
                - (2.0 * (jf.j2s - cache.ys * h2s / 2.0) + jf.j2c - cache.ys * h2c / 2.0) / 3.0,
            2.0 * (jt.j1s - cache.ys * h1s / 2.0) + jt.j1c - cache.ys * h1c / 2.0
                - (2.0 * (jt.j2s - cache.ys * h2s / 2.0) + jt.j2c - cache.ys * h2c / 2.0) / 3.0,
            jt.j2s + jf.j2s - cache.ys * h2s,
            jt.j2c + jf.j2c - cache.ys * h2c,
            jt.j3  + jf.j3  - cache.ys * h3,
            jt.j4  + jf.j4  - cache.ys * h4,
            val6,
            val7,
            beta * val7,
            val9,
            jt.j7  + jf.j7  - cache.ys * h7,
            val11,
            val12,
            s8v,
            s9v
        }};
    };

    auto integrate_bin = [&] (double q2_l, double q2_u) -> std::array<double, 15> {
        std::array<double, 15> acc {};
        const double center = 0.5 * (q2_l + q2_u);
        const double half_width = 0.5 * (q2_u - q2_l);

        for (size_t i = 0; i < GL24_X.size(); ++i) {
            const double q2 = center + half_width * GL24_X[i];
            const auto vals = eval_integrands(q2);
            for (size_t k = 0; k < acc.size(); ++k) {
                acc[k] += GL24_W[i] * vals[k];
            }
        }

        for (double& v : acc) {
            v *= half_width;
        }

        return acc;
    };

    const double endpoint_eps = 1e-5;
    const double low_q2_eps = 1e-7;

    const auto& bins = this->bins.value();
    const size_t nbins = bins.size();
    clear_and_reserve(cache.f_J_i_binned, nbins);
    cache.bin_widths.clear();
    cache.bin_widths.reserve(nbins);

    for (const auto& [q2_l_raw, q2_u_raw] : bins) {
        const double q2_l = std::max(q2_l_raw, cache.q2_min + low_q2_eps);
        const double q2_u = std::min(q2_u_raw, cache.q2_max - endpoint_eps);

        if (!(q2_l < q2_u)) {
            LOG_WARN(
                "Skipping invalid BsPhi bin [",
                q2_l_raw,
                ",",
                q2_u_raw,
                "] clipped to [",
                q2_l,
                ",",
                q2_u,
                "]"
            );

            for (auto& v : cache.f_J_i_binned) {
                v.emplace_back(std::numeric_limits<double>::quiet_NaN());
            }
            cache.bin_widths.emplace_back(std::numeric_limits<double>::quiet_NaN());
            continue;
        }

        const double width = q2_u - q2_l;
        const auto integ = integrate_bin(q2_l, q2_u);

        cache.bin_widths.emplace_back(width);

        for (size_t k = 0; k < cache.f_J_i_binned.size(); ++k) {
            cache.f_J_i_binned[k].emplace_back(integ[k]);
        }
    }
}

std::vector<ObservableValue> BsPhiDecay::dBR_dq2_binned(bool bar, Observables id, bool br) {
    std::vector<ObservableValue> out;
    size_t idx = bar ? 1 : 0;
    const auto& bins = this->bins.value();
    const double br_factor = br ? cache.life_Bs : 1.0;

    for (size_t i = 0; i < bins.size(); i++) {
        const double requested_width = bins[i].second - bins[i].first;
        const double width =
            (i < cache.bin_widths.size() && std::isfinite(cache.bin_widths[i]) && cache.bin_widths[i] > 0.0)
            ? cache.bin_widths[i]
            : requested_width;

        const double integrated_rate =
            0.75 * cache.f_J_i_binned[idx][i] / (1.0 - cache.ys * cache.ys);

        const double res =
            (std::isfinite(width) && width > 0.0)
            ? integrated_rate * br_factor / width
            : std::numeric_limits<double>::quiet_NaN();

        out.emplace_back(ObservableMapper::to_id(id), res, bins[i]);
    }

    return out;
}

double BsPhiDecay::dG_dq2_avg_bin(size_t bin) {
    return 0.75 * (cache.f_J_i_binned[0][bin] + cache.f_J_i_binned[1][bin]);
}

static double safe_ratio(double num, double den, const char* name) {
    if (!std::isfinite(num) || !std::isfinite(den) || std::abs(den) < 1e-30) {
        LOG_WARN(
            "Non-finite or singular denominator in",
            name,
            ": num =",
            num,
            "den =",
            den
        );
        return std::numeric_limits<double>::quiet_NaN();
    }

    return num / den;
}

std::vector<ObservableValue> BsPhiDecay::F_L(Observables id) {
    std::vector<ObservableValue> out;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        const double res = safe_ratio(
            -cache.f_J_i_binned[3][i],
            dG_dq2_avg_bin(i),
            "BsPhi F_L"
        );

        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }

    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_T_2(Observables id) {
    std::vector<ObservableValue> out;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        const double res = safe_ratio(
            0.5 * cache.f_J_i_binned[4][i],
            cache.f_J_i_binned[2][i],
            "BsPhi A_T_2"
        );

        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }

    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_T_Re_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = 0.25 * cache.f_J_i_binned[8][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_T_Im_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[12][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_4(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = cache.f_J_i_binned[5][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_6(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = -0.5 * cache.f_J_i_binned[10][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::S_i(int i, Observables id) {
    if (!(i == 2 || i == 3 || i == 4 || i == 7)) LOG_ERROR("Value Error", "S_i(Bs > phi ll) is not defined for i =", i);

    std::map<size_t, size_t> J_idx = {{2, 2}, {3, 4}, {4, 5}, {7, 10}};

    std::vector<ObservableValue> out;
    for (size_t j = 0; j < this->bins.value().size(); j++) {
        double res = cache.f_J_i_binned[J_idx[i]][j] / dG_dq2_avg_bin(j);
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[j]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_i(int i, Observables id) {
    if (!(i == 5 || i == 6 || i == 8 || i == 9)) LOG_ERROR("Value Error", "A_i(Bs > phi ll) is not defined for i =", i);

    std::map<size_t, size_t> J_idx = {{5, 6}, {6, 9}, {8, 11}, {9, 12}};
    double sign = i == 6 ? -1 : 1;

    std::vector<ObservableValue> out;
    for (size_t j = 0; j < this->bins.value().size(); j++) {
        double res = sign * cache.f_J_i_binned[J_idx[i]][j] / dG_dq2_avg_bin(j);
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[j]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::A_FB_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = -0.375 * (2 * cache.f_J_i_binned[7][i] + cache.f_J_i_binned[9][i]) / dG_dq2_avg_bin(i); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::P_2_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = 0.125 * cache.f_J_i_binned[7][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::P_3_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = -0.25 * cache.f_J_i_binned[12][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_5_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[6][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Pp_8_CPV(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = -cache.f_J_i_binned[11][i] / std::sqrt(std::abs(cache.f_J_i_binned[2][i] * cache.f_J_i_binned[3][i])); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Q_8_m(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = cache.f_J_i_binned[13][i] / std::sqrt(std::abs(2 * cache.f_J_i_binned[3][i] * (cache.f_J_i_binned[2][i] - cache.f_J_i_binned[4][i]))); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Q_8_p(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = cache.f_J_i_binned[13][i] / std::sqrt(std::abs(2 * cache.f_J_i_binned[3][i] * (cache.f_J_i_binned[2][i] + cache.f_J_i_binned[4][i]))); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Q_9(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = 0.5 * cache.f_J_i_binned[14][i] / cache.f_J_i_binned[2][i]; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BsPhiDecay::Rm1_BsPhi(Observables id) {
    std::vector<ObservableValue> out;
    std::vector<double> Gamma_mu;
    std::vector<double> Gamma_e;

    set_cfg_flags(BsPhiConfig::Lepton::MU);

    for (size_t i = 0; i < this->bins.value().size(); i++)
        Gamma_mu.emplace_back(dG_dq2_avg_bin(i));

    set_cfg_flags(BsPhiConfig::Lepton::E);

    for (size_t i = 0; i < this->bins.value().size(); i++)
        Gamma_e.emplace_back(dG_dq2_avg_bin(i));

    for (size_t i = 0; i < this->bins.value().size(); i++)
        out.emplace_back(ObservableMapper::to_id(id), Gamma_mu[i] / Gamma_e[i] - 1, this->bins.value()[i]);
    
    return out;
}

std::vector<ObservableValue> BsPhiDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::DGAMMA_DQ2_BS__PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return dBR_dq2_binned(false, obs, false);
    case Observables::DBR_DQ2_BS__PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return dBR_dq2_binned(false, obs);
    case Observables::F_L_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return F_L(obs);
    case Observables::A_T_2_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return A_T_2(obs);
    case Observables::A_T_RE_CPV_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return A_T_Re_CPV(obs);
    case Observables::A_T_IM_CPV_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return A_T_Im_CPV(obs);
    case Observables::P_PRIME_4_BS_PHI_E_E:      
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return Pp_4(obs);
    case Observables::P_PRIME_6_BS_PHI_E_E:      
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return Pp_6(obs);
    case Observables::S_2S_BS_PHI_E_E:      
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return S_i(2, obs);
    case Observables::S_3_BS_PHI_E_E:      
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return S_i(3, obs);
    case Observables::S_4_BS_PHI_E_E:      
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return S_i(4, obs);
    case Observables::S_7_BS_PHI_E_E:      
        set_cfg_flags(BsPhiConfig::Lepton::E);
        return S_i(7, obs);
    case Observables::A_5_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return A_i(5, obs);
    case Observables::A_6C_BS_PHI_E_E:     
        set_cfg_flags(BsPhiConfig::Lepton::E); 
        return A_i(6, obs);
    case Observables::A_8_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return A_i(8, obs);
    case Observables::A_9_BS_PHI_E_E:    
        set_cfg_flags(BsPhiConfig::Lepton::E);  
        return A_i(9, obs);
    case Observables::A_FB_CPV_BS__PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return A_FB_CPV(obs);
    case Observables::P_2_CPV_BS_PHI_E_E:    
        set_cfg_flags(BsPhiConfig::Lepton::E);  
        return P_2_CPV(obs);
    case Observables::P_3_CPV_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return P_3_CPV(obs);
    case Observables::P_PRIME_5_CPV_BS_PHI_E_E:   
        set_cfg_flags(BsPhiConfig::Lepton::E);   
        return Pp_5_CPV(obs);
    case Observables::P_PRIME_8_CPV_BS_PHI_E_E:    
        set_cfg_flags(BsPhiConfig::Lepton::E);  
        return Pp_8_CPV(obs);
    case Observables::Q_8M_BS_PHI_E_E:     
        set_cfg_flags(BsPhiConfig::Lepton::E); 
        return Q_8_m(obs);
    case Observables::Q_8P_BS_PHI_E_E:     
        set_cfg_flags(BsPhiConfig::Lepton::E); 
        return Q_8_p(obs);
    case Observables::Q_9_BS_PHI_E_E:    
        set_cfg_flags(BsPhiConfig::Lepton::E);  
        return Q_9(obs);
    case Observables::DGAMMA_DQ2_BS__PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return dBR_dq2_binned(false, obs, false);
    case Observables::DBR_DQ2_BS__PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return dBR_dq2_binned(false, obs);
    case Observables::F_L_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return F_L(obs);
    case Observables::A_T_2_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return A_T_2(obs);
    case Observables::A_T_RE_CPV_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return A_T_Re_CPV(obs);
    case Observables::A_T_IM_CPV_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return A_T_Im_CPV(obs);
    case Observables::P_PRIME_4_BS_PHI_MU_MU:      
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return Pp_4(obs);
    case Observables::P_PRIME_6_BS_PHI_MU_MU:      
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return Pp_6(obs);
    case Observables::S_2S_BS_PHI_MU_MU:      
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return S_i(2, obs);
    case Observables::S_3_BS_PHI_MU_MU:      
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return S_i(3, obs);
    case Observables::S_4_BS_PHI_MU_MU:      
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return S_i(4, obs);
    case Observables::S_7_BS_PHI_MU_MU:      
        set_cfg_flags(BsPhiConfig::Lepton::MU);
        return S_i(7, obs);
    case Observables::A_5_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return A_i(5, obs);
    case Observables::A_6C_BS_PHI_MU_MU:     
        set_cfg_flags(BsPhiConfig::Lepton::MU); 
        return A_i(6, obs);
    case Observables::A_8_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return A_i(8, obs);
    case Observables::A_9_BS_PHI_MU_MU:    
        set_cfg_flags(BsPhiConfig::Lepton::MU);  
        return A_i(9, obs);
    case Observables::A_FB_CPV_BS__PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return A_FB_CPV(obs);
    case Observables::P_2_CPV_BS_PHI_MU_MU:    
        set_cfg_flags(BsPhiConfig::Lepton::MU);  
        return P_2_CPV(obs);
    case Observables::P_3_CPV_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return P_3_CPV(obs);
    case Observables::P_PRIME_5_CPV_BS_PHI_MU_MU:   
        set_cfg_flags(BsPhiConfig::Lepton::MU);   
        return Pp_5_CPV(obs);
    case Observables::P_PRIME_8_CPV_BS_PHI_MU_MU:    
        set_cfg_flags(BsPhiConfig::Lepton::MU);  
        return Pp_8_CPV(obs);
    case Observables::Q_8M_BS_PHI_MU_MU:     
        set_cfg_flags(BsPhiConfig::Lepton::MU); 
        return Q_8_m(obs);
    case Observables::Q_8P_BS_PHI_MU_MU:     
        set_cfg_flags(BsPhiConfig::Lepton::MU); 
        return Q_8_p(obs);
    case Observables::Q_9_BS_PHI_MU_MU:    
        set_cfg_flags(BsPhiConfig::Lepton::MU);  
        return Q_9(obs);
    case Observables::DGAMMA_DQ2_BS__PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return dBR_dq2_binned(false, obs, false);
    case Observables::DBR_DQ2_BS__PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return dBR_dq2_binned(false, obs);
    case Observables::F_L_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return F_L(obs);
    case Observables::A_T_2_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return A_T_2(obs);
    case Observables::A_T_RE_CPV_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return A_T_Re_CPV(obs);
    case Observables::A_T_IM_CPV_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return A_T_Im_CPV(obs);
    case Observables::P_PRIME_4_BS_PHI_TAU_TAU:      
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return Pp_4(obs);
    case Observables::P_PRIME_6_BS_PHI_TAU_TAU:      
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return Pp_6(obs);
    case Observables::S_2S_BS_PHI_TAU_TAU:      
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return S_i(2, obs);
    case Observables::S_3_BS_PHI_TAU_TAU:      
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return S_i(3, obs);
    case Observables::S_4_BS_PHI_TAU_TAU:      
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return S_i(4, obs);
    case Observables::S_7_BS_PHI_TAU_TAU:      
        set_cfg_flags(BsPhiConfig::Lepton::TAU);
        return S_i(7, obs);
    case Observables::A_5_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return A_i(5, obs);
    case Observables::A_6C_BS_PHI_TAU_TAU:     
        set_cfg_flags(BsPhiConfig::Lepton::TAU); 
        return A_i(6, obs);
    case Observables::A_8_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return A_i(8, obs);
    case Observables::A_9_BS_PHI_TAU_TAU:    
        set_cfg_flags(BsPhiConfig::Lepton::TAU);  
        return A_i(9, obs);
    case Observables::A_FB_CPV_BS__PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return A_FB_CPV(obs);
    case Observables::P_2_CPV_BS_PHI_TAU_TAU:    
        set_cfg_flags(BsPhiConfig::Lepton::TAU);  
        return P_2_CPV(obs);
    case Observables::P_3_CPV_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return P_3_CPV(obs);
    case Observables::P_PRIME_5_CPV_BS_PHI_TAU_TAU:   
        set_cfg_flags(BsPhiConfig::Lepton::TAU);   
        return Pp_5_CPV(obs);
    case Observables::P_PRIME_8_CPV_BS_PHI_TAU_TAU:    
        set_cfg_flags(BsPhiConfig::Lepton::TAU);  
        return Pp_8_CPV(obs);
    case Observables::Q_8M_BS_PHI_TAU_TAU:     
        set_cfg_flags(BsPhiConfig::Lepton::TAU); 
        return Q_8_m(obs);
    case Observables::Q_8P_BS_PHI_TAU_TAU:     
        set_cfg_flags(BsPhiConfig::Lepton::TAU); 
        return Q_8_p(obs);
    case Observables::Q_9_BS_PHI_TAU_TAU:    
        set_cfg_flags(BsPhiConfig::Lepton::TAU);  
        return Q_9(obs);
    case Observables::R_1_BS__PHI_L_L:
        return Rm1_BsPhi(obs);
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> BsPhiDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}

void BsPhiDecay::set_n_threads(size_t n_threads) {
    unsigned int available_threads = std::thread::hardware_concurrency();

    if (available_threads == 0) {
        available_threads = 1;
    }

    if (n_threads == 0) {
        this->cfg.n_threads = available_threads;
        return;
    }

    if (n_threads > available_threads) {
        LOG_WARN(
            "Requested", n_threads,
            "threads, but only", available_threads,
            "are available. Using", available_threads,
            "threads instead."
        );

        this->cfg.n_threads = available_threads;
        return;
    }

    this->cfg.n_threads = std::max<size_t>(1, n_threads);
}
