
#include "BKsllDecay.h"
#include "wcoef_ids.hpp"

#include <algorithm>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

using Charge = BKstarllConfig::B_Charge;

static double finite_or_nan(double x) {
    return std::isfinite(x) ? x : std::numeric_limits<double>::quiet_NaN();
}

static double safe_div(double num, double den, const char* name) {
    if (!std::isfinite(num) || !std::isfinite(den) || std::abs(den) < 1e-30) {
        LOG_WARN(name, "non-finite/singular ratio: num =", num, "den =", den);
        return std::numeric_limits<double>::quiet_NaN();
    }

    return num / den;
}

void BKstarllDecay::load_params() {
    fill_wilson_cache();

    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.G_F = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    cache.m_b_mb = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 5}, DataType::VALUE);
    cache.m_s = (*p)(ParamId{ParameterType::SM, "MASS", 3}, DataType::VALUE);
    cache.mu_b = (*p)(ParamId{ParameterType::WILSON, "B_SCALE", 1}, DataType::VALUE);
    cache.alpha_s_mu_b = (*iobs_qcdp)(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    // cache.m_c_mu_b = (*ports.iobs_qcdp)(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_c_mu_b = (*p)(ParamId{ParameterType::SM, "MASS", 4}, DataType::VALUE); // ASK : To match SI, should probably be m_c(mu_b) instead
    cache.m_b_mu_b = (*iobs_qcdp)(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    double mu_f = sqrt(cache.mu_b * (*p)(ParamId{ParameterType::DECAY, "B_Ks", 14}, DataType::VALUE));
    cache.m_b_PS = (*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE) - 4 * (*iobs_qcdp)(AlphasConfig((*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE), MassType::POLE, MassType::POLE)) * mu_f / (3 * PI);
    cache.L_b = std::log(cache.mu_b / cache.m_b_PS);
    cache.Delta_M = -6. * cache.L_b - 4. * (1 - mu_f / cache.m_b_PS);
    cache.lambda_hat_u = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE) 
                            / (std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE));
    cache.kappa = 1. - 2. * cache.alpha_s_mu_b / (3. * PI) * std::log(cache.mu_b / cache.m_b_mb);
    cache.q2_low = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {15, 1}}, DataType::VALUE);
    cache.q2_high = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {15, 2}}, DataType::VALUE);

    for (size_t i = 0; i < 6; i++) {
        cache.a_k_low[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {18, 1, i + 1}}, DataType::VALUE);
        cache.b_k_low[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {18, 2, i + 1}}, DataType::VALUE);
        cache.phi_k_low[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {18, 4, i + 1}}, DataType::VALUE);
        cache.theta_k_low[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {18, 5, i + 1}}, DataType::VALUE);
    }

    for (size_t i = 0; i < 8; i++) {
        cache.a_k_high[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {18, 3, i + 1}}, DataType::VALUE);
        cache.phi_k_high[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {18, 6, i + 1}}, DataType::VALUE);
    }

    load_cfg_dependent_params();
}

// void BKstarllDecay::fill_wilson_cache() {
//     auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
//     auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
//     auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);
//     // All BPrime coefficients are cached because CPQ1/CPQ2 are lepton-specific.

//     for (auto p : b_wilsons) cache.C.emplace(p); 
//     for (auto p : bq_wilsons) cache.C.emplace(p);
//     for (auto id : bp_cached) cache.C.emplace(std::pair{id, bp_wilsons.at(id)});
// }

// void BKstarllDecay::fill_wilson_cache() {
//     auto b_wilsons  = w_proxy->getAFR(WGroup::B, this->w_config.order);
//     auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
//     auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);

//     static std::size_t dbg = 0;
//     if (dbg < 20) {
//         std::cout << "[DECDBG] fill_wilson_cache call=" << dbg << "\n";
//         std::cout << "[DECDBG] B:C9   = " << b_wilsons.at(WCoef::C9) << "\n";
//         std::cout << "[DECDBG] B:C10  = " << b_wilsons.at(WCoef::C10) << "\n";
//         std::cout << "[DECDBG] BP:CP9 = " << bp_wilsons.at(WCoef::CP9) << "\n";
//         std::cout << "[DECDBG] BP:CP10= " << bp_wilsons.at(WCoef::CP10) << "\n";
//     }
//     ++dbg;

//     // All BPrime coefficients are cached because CPQ1/CPQ2 are lepton-specific.

//     for (auto p : b_wilsons) cache.C.emplace(p);
//     for (auto p : bq_wilsons) cache.C.emplace(p);
//     for (auto id : bp_cached) cache.C.emplace(std::pair{id, bp_wilsons.at(id)});
// }

void BKstarllDecay::fill_wilson_cache() {
    auto b_wilsons  = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);

    // static std::size_t dbg = 0;
    // if (dbg < 20) {
    //     std::cout << "[DECDBG] fill_wilson_cache call=" << dbg << "\n";
    //     std::cout << "[DECDBG] source B:C9   = " << b_wilsons.at(WCoef::C9) << "\n";
    //     std::cout << "[DECDBG] source B:C10  = " << b_wilsons.at(WCoef::C10) << "\n";
    //     std::cout << "[DECDBG] source BP:CP9 = " << bp_wilsons.at(WCoef::CP9) << "\n";
    //     std::cout << "[DECDBG] source BP:CP10= " << bp_wilsons.at(WCoef::CP10) << "\n";
    // }

    // CRUCIAL: on repart d'un cache vide, sinon emplace garde les vieilles valeurs
    cache.C.clear();

    for (const auto& [coef, val] : b_wilsons) {
        cache.C[coef] = val;
    }
    for (const auto& [coef, val] : bq_wilsons) {
        cache.C[coef] = val;
    }

    for (const auto& [coef, val] : bp_wilsons) {
        cache.C[coef] = val;
    }

    // if (dbg < 20) {
    //     std::cout << "[DECDBG] cache  B:C9   = " << cache.C.at(WCoef::C9) << "\n";
    //     std::cout << "[DECDBG] cache  B:C10  = " << cache.C.at(WCoef::C10) << "\n";
    //     std::cout << "[DECDBG] cache  BP:CP9 = " << cache.C.at(WCoef::CP9) << "\n";
    //     std::cout << "[DECDBG] cache  BP:CP10= " << cache.C.at(WCoef::CP10) << "\n";
    // }
    // ++dbg;
}

void BKstarllDecay::load_cfg_dependent_params() {
    const int B_id = cfg.charge == Charge::B_0 ? 511 : 521;
    const int V_id = cfg.charge == Charge::B_0 ? 313 : 323;

    cache.ff_calculator = BVFFCalculator(
        B_id,
        V_id,
        p,
        cfg.ff_src
    );

    cache.qcdf_calculator = BVQCDfCalculator(
        B_id,
        V_id,
        cache.mu_b,
        cache.C,
        std::make_shared<BVFFCalculator>(cache.ff_calculator),
        cfg.ff_type,
        p,
        iobs_qcdp
    );

    cache.m_l = (*p)(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen}, DataType::VALUE);
    cache.m_B = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 511 : 521}, DataType::VALUE);
    cache.life_B = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", cfg.charge == Charge::B_0 ? 511 : 521}, DataType::VALUE) / HBAR;
    cache.m_Ks = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 313 : 323}, DataType::VALUE);
    cache.N_0 = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * cache.G_F * cache.alpha_em / (std::sqrt(3072. * std::pow(PI, 5) * std::pow(cache.m_B, 3)));
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    cache.q2_max = std::pow(cache.m_B - cache.m_Ks, 2);

    cache.q2_lookup_min = std::max(cache.q2_min, 1e-4);
    if (cfg.ff_type == B_FF_Type::SOFT || cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) {
        auto requested_threads = cfg.n_threads;
        if (requested_threads == 0u) {
            requested_threads = std::thread::hardware_concurrency();
        }
        if (requested_threads == 0u) {
            requested_threads = 1u;
        }

        const size_t npts = BKstarllCache::LOOKUP_SIZE;
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
    }

    if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BCvDV) {
        double m_D0 = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 421}, DataType::VALUE);
        double m_Jpsi = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 443}, DataType::VALUE);
        double m_psi_2S = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 100443}, DataType::VALUE);
        cache.tp_nf = 4. * std::pow(m_D0, 2);
        cache.t0_nf = cache.tp_nf - std::sqrt(cache.tp_nf * (cache.tp_nf - std::pow(m_psi_2S, 2)));
        cache.z0_nf = (sqrt(cache.tp_nf) - sqrt(cache.tp_nf - cache.t0_nf)) / (sqrt(cache.tp_nf) + sqrt(cache.tp_nf - cache.t0_nf));
        cache.z_Jpsi_nf = (sqrt(cache.tp_nf - std::pow(m_Jpsi, 2)) - sqrt(cache.tp_nf - cache.t0_nf)) / (sqrt(cache.tp_nf - std::pow(m_Jpsi, 2)) + sqrt(cache.tp_nf - cache.t0_nf));
        cache.z_psi2S_nf = (sqrt(cache.tp_nf - std::pow(m_psi_2S, 2)) - sqrt(cache.tp_nf - cache.t0_nf)) / (sqrt(cache.tp_nf - std::pow(m_psi_2S, 2)) + sqrt(cache.tp_nf - cache.t0_nf));

        for (size_t i = 0; i < 3; i++) {
            cache.alpha_perp[i] = complex_t {
                (*p)(ParamId{ParameterType::DECAY, "B_Ks", {19, 1, 1, i}}, DataType::VALUE),
                (*p)(ParamId{ParameterType::DECAY, "B_Ks", {19, 2, 1, i}}, DataType::VALUE)
            };

            cache.alpha_par[i] = complex_t {
                (*p)(ParamId{ParameterType::DECAY, "B_Ks", {19, 1, 2, i}}, DataType::VALUE),
                (*p)(ParamId{ParameterType::DECAY, "B_Ks", {19, 2, 2, i}}, DataType::VALUE)
            };
        }
        
        for (size_t i = 0; i < 2; i++) {
            cache.alpha_0[i] = complex_t {
                (*p)(ParamId{ParameterType::DECAY, "B_Ks", {19, 1, 3, i}}, DataType::VALUE),
                (*p)(ParamId{ParameterType::DECAY, "B_Ks", {19, 2, 3, i}}, DataType::VALUE)
            };
        }
    }

    if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::KMPW) {
        cache.q2_bar = 1.0;
        cache.q2_Jpsi = std::pow(std::real((*p)(ParamId{ParameterType::FLAVOR, "FMASS", 443}, DataType::VALUE)), 2);

        for (size_t i = 0; i < 3; i++) {
            cache.DeltaC9_M_qbar[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {20, i}}, DataType::VALUE);
            cache.r1_M[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {21, 1, i}}, DataType::VALUE);
            cache.r2_M[i] = (*p)(ParamId{ParameterType::DECAY, "B_Ks", {21, 2, i}}, DataType::VALUE);
        }
    }

    // double q2 = 14.0;

    // printf("alpha_em = %.4e\n", cache.alpha_em);
    // printf("kappa = %.4e\n", cache.kappa);
    // printf("N_0 = %.4e + %.4e i\n", std::real(cache.N_0 * std::sqrt(q2 * std::sqrt(lambda(q2)) * beta_l(q2))), std::imag(cache.N_0 * std::sqrt(q2 * std::sqrt(lambda(q2)) * beta_l(q2))));

    // double sh = q2 / std::pow(cache.m_b_PS, 2);

    // printf("mb_PS = %.4e\n", cache.m_b_PS);
    // printf("A = %.4e + %.4e i\n", real(BV::A_Seidel(sh, cache.L_b)), imag(BV::A_Seidel(sh, cache.L_b)));
    // printf("B = %.4e + %.4e i\n", real(BV::B_Seidel(sh, cache.L_b)), imag(BV::B_Seidel(sh, cache.L_b)));
    // printf("C = %.4e + %.4e i\n", real(BV::C_Seidel(q2, cache.mu_b)), imag(BV::C_Seidel(q2, cache.mu_b)));
    // printf("alpha_s_mu_b = %.4e\n", cache.alpha_s_mu_b);
    // printf("m_b_mu_b = %.4e\n", cache.m_b_mu_b);
    // printf("kappa = %.4e\n", cache.kappa);

    // LOG_INFO("Anchor");
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
    // exit(0);
    
    test_J();

    compute_binned_J_i();
}

void BKstarllDecay::set_lepton_gen_and_charge(BKstarllConfig::Lepton gen, BKstarllConfig::B_Charge charge) {
    bool changed = cfg.gen != gen || cfg.charge != charge;
    if (changed) {
        cfg.gen = gen;  
        cfg.charge = charge;
        load_cfg_dependent_params();
    }
}

complex_t BKstarllDecay::T_perp_p_cached(double q2, bool bar) {
    const double x = std::clamp(q2, cache.q2_lookup_min, cache.q2_high);
    return lerp(
        x,
        bar ? cache.T_perp_p_bar_lookup : cache.T_perp_p_lookup,
        cache.q2_lookup_min,
        cache.q2_high
    );
}

complex_t BKstarllDecay::T_perp_m_cached(double q2, bool bar) {
    const double x = std::clamp(q2, cache.q2_lookup_min, cache.q2_high);
    return lerp(
        x,
        bar ? cache.T_perp_m_bar_lookup : cache.T_perp_m_lookup,
        cache.q2_lookup_min,
        cache.q2_high
    );
}

complex_t BKstarllDecay::T_par_m_cached(double q2, bool bar) {
    const double x = std::clamp(q2, cache.q2_lookup_min, cache.q2_high);
    return lerp(
        x,
        bar ? cache.T_par_m_bar_lookup : cache.T_par_m_lookup,
        cache.q2_lookup_min,
        cache.q2_high
    );
}

double BKstarllDecay::beta_l(double q2) {
    const double x = 1.0 - 4.0 * cache.m_l * cache.m_l / q2;
    return std::sqrt(std::max(0.0, x));
}

double BKstarllDecay::lambda(double q2) {
    const double mB2 = cache.m_B * cache.m_B;
    const double mK2 = cache.m_Ks * cache.m_Ks;

    const double lam =
        mB2 * mB2
        + mK2 * mK2
        + q2 * q2
        - 2.0 * (mB2 * mK2 + (mB2 + mK2) * q2);

    return std::max(0.0, lam);
}

complex_t BKstarllDecay::N(double q2, bool bar) {
    const complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0;
    const double b = beta_l(q2);
    const double lam = lambda(q2);

    return N0 * std::sqrt(std::max(0.0, q2 * b * std::sqrt(lam)));
}

complex_t BKstarllDecay::delta_A_perp_QCDf(double q2, double sign, bool bar) {
    complex_t guesstimate_err = 1.0;
    complex_t delta_A = 0.0;

    if (!fpeq(std::abs(cache.h_p_fit[0]), 0.0)) {
        complex_t h_p = cache.h_p_fit[0] + q2 * (cache.h_p_fit[1] + q2 * cache.h_p_fit[2]);
        complex_t h_m = cache.h_m_fit[0] + q2 * (cache.h_m_fit[1] + q2 * cache.h_m_fit[2]);
        delta_A = 16.0 * PI2 * RT2 * N(q2, bar) * std::pow(cache.m_B, 3) * (h_p - h_m) / q2;
    } else {
        size_t id = size_t (0.5 * (1 + sign));
        guesstimate_err = 1.0 + cache.a_k_low[id] * std::exp(I * cache.phi_k_low[id]) + cache.b_k_low[id] * std::exp(I * cache.theta_k_low[id]) * q2 / 6.0;
    }

    // LOG_INFO("Guesstimate error in delta_A_perp_QCDf = ", guesstimate_err);

    return 2 * RT2 * (cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI)) * N(q2, bar) * std::sqrt(lambda(q2)) / q2 * T_perp_p_cached(q2, bar) * guesstimate_err + delta_A;
}

complex_t BKstarllDecay::delta_A_perp_vD(double q2, bool bar) {
    double F_perp = std::sqrt(2. * lambda(q2)) / (cache.m_B * (cache.m_B + cache.m_Ks)) * cache.ff_calculator.get(BV_FF::V, q2);
    complex_t z_q2 = cache.ff_calculator.z(q2, cache.tp_nf, cache.t0_nf);
    complex_t P_H_perp = cache.alpha_perp[0] + z_q2 * (cache.alpha_perp[1] + z_q2 * cache.alpha_perp[2]);
    complex_t H_perp = (1. - z_q2 * std::conj(cache.z_Jpsi_nf)) / (z_q2 - cache.z_Jpsi_nf) * (1. - z_q2 * std::conj(cache.z_psi2S_nf)) / (z_q2 - cache.z_psi2S_nf) * P_H_perp * F_perp;
    return -32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 3) * H_perp / q2;
}

complex_t BKstarllDecay::delta_A_perp_K(double q2, bool bar) {
    double DeltaC9_M1 = (cache.r1_M[0] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[0] * cache.q2_bar / q2) / (1 + cache.r2_M[0] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    return N(q2, bar) * RT2 * std::sqrt(lambda(q2)) * DeltaC9_M1 * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_B + cache.m_Ks);
}

complex_t BKstarllDecay::delta_A_perp(double q2, double sign, bool bar) {
    switch(cfg.power_corr_impl) {
        case BKstarllConfig::Power_Corrections_Impl::BFS:
            return delta_A_perp_QCDf(q2, sign, bar);
        case BKstarllConfig::Power_Corrections_Impl::BCvDV:
            return delta_A_perp_vD(q2, bar);
        case BKstarllConfig::Power_Corrections_Impl::KMPW:
            return delta_A_perp_K(q2, bar);
        default:
            return 0.0;
    }
}

complex_t BKstarllDecay::delta_A_par_QCDf(double q2, double sign, bool bar) {
    complex_t guesstimate_err = 1.0;
    complex_t delta_A = 0.0;

    if (!fpeq(std::abs(cache.h_p_fit[0]), 0.0)) {
        complex_t h_p = cache.h_p_fit[0] + q2 * (cache.h_p_fit[1] + q2 * cache.h_p_fit[2]);
        complex_t h_m = cache.h_m_fit[0] + q2 * (cache.h_m_fit[1] + q2 * cache.h_m_fit[2]);
        delta_A = 16.0 * PI2 * RT2 * N(q2, bar) * std::pow(cache.m_B, 3) * (h_p + h_m) / q2;
    } else {
        size_t id = 2 + size_t (0.5 * (1 + sign));
        guesstimate_err = guesstimate_err = 1.0 + cache.a_k_low[id] * std::exp(I * cache.phi_k_low[id]) + cache.b_k_low[id] * std::exp(I * cache.theta_k_low[id]) * q2 / 6.0;
    }

    return -4 * RT2 * (cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI)) * N(q2, bar) * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks) * cache.ff_calculator.E(q2) / (q2 * cache.m_B) * T_perp_m_cached(q2, bar) * guesstimate_err + delta_A;
}

complex_t BKstarllDecay::delta_A_par_vD(double q2, bool bar) {
    double F_par = RT2 * (cache.m_B + cache.m_Ks) / cache.m_B * cache.ff_calculator.get(BV_FF::A1, q2);
    complex_t z_q2 = cache.ff_calculator.z(q2, cache.tp_nf, cache.t0_nf);
    complex_t P_H_par = cache.alpha_par[0] + z_q2 * (cache.alpha_par[1] + z_q2 * cache.alpha_par[2]);
    complex_t H_par = (1. - z_q2 * std::conj(cache.z_Jpsi_nf)) / (z_q2 - cache.z_Jpsi_nf) * (1. - z_q2 * std::conj(cache.z_psi2S_nf)) / (z_q2 - cache.z_psi2S_nf) * P_H_par * F_par;
    return 32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 3) * H_par / q2;
}

complex_t BKstarllDecay::delta_A_par_K(double q2, bool bar) {
    double DeltaC9_M2 = (cache.r1_M[1] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[1] * cache.q2_bar / q2) / (1 + cache.r2_M[1] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    return N(q2, bar) / RT2 * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks) * DeltaC9_M2 * cache.ff_calculator.get(BV_FF::A1, q2) / (cache.m_B - cache.m_Ks);
}

complex_t BKstarllDecay::A_perp_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_B + cache.m_Ks);
        F_T = T_perp_p_cached(q2, bar);
        size_t id = size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.a_k_low[id] + cache.b_k_low[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BV_FF::T1, q2);
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) w += cache.qcdf_calculator.Y(q2);
        F = w * cache.ff_calculator.get(BV_FF::V, q2) / (cache.m_B + cache.m_Ks);
        delta_A = delta_A_perp(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("----- A_perp -----\n");
    // printf("pref = %.4e + %.4e i\n", real(N(q2, bar) * std::sqrt(2 * lambda(q2))), imag(N(q2, bar) * std::sqrt(2 * lambda(q2))));
    // printf("F_V = %.4e + %.4e i\n", real(F), imag(F));
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));
    // printf("delta_A = %.4e + %.4e i\n", real(delta_A), imag(delta_A));

    return (N(q2, bar) * std::sqrt(2 * lambda(q2)) * (F + 2. * m_b_local * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BKstarllDecay::A_par_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * cache.ff_calculator.get(BV_FF::XI_PERP, q2) * 2. * cache.ff_calculator.E(q2) / (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks);
        F_T = 2. * cache.ff_calculator.E(q2) * T_perp_m_cached(q2, bar) / cache.m_B;
        size_t id = 2 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.a_k_low[id] + cache.b_k_low[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * cache.ff_calculator.get(BV_FF::T2, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) w += cache.qcdf_calculator.Y(q2);
        F = w * cache.ff_calculator.get(BV_FF::A1, q2) / (cache.m_B - cache.m_Ks);
        delta_A = delta_A_par(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("----- A_par -----\n");
    // printf("pref = %.4e + %.4e i\n", real(-N(q2, bar) * std::sqrt(2.) * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks)), imag(-N(q2, bar) * std::sqrt(2.) * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks)));
    // printf("F_V = %.4e + %.4e i\n", real(F), imag(F));
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));
    // printf("delta_A = %.4e + %.4e i\n", real(delta_A), imag(delta_A));

    return (-N(q2, bar) * std::sqrt(2.) * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks) * (F + 2. * m_b_local * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BKstarllDecay::A_0_low(double q2, double sign, bool bar) {
    double mB2 = cache.m_B * cache.m_B;
    double mK2 = cache.m_Ks * cache.m_Ks;
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};
    double m_b_local;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * (
            (2 * cache.ff_calculator.E(q2) * (mB2 - mK2 - q2) - lambda(q2) * cache.m_B / (mB2 - mK2)) * cache.ff_calculator.get(BV_FF::XI_PERP, q2) 
          + (lambda(q2) * cache.m_B / (mB2 - mK2)) * cache.ff_calculator.get(BV_FF::XI_PAR, q2)
        );
        F_T = 2. * cache.ff_calculator.E(q2) * (mB2 + 3. * mK2 - q2) / cache.m_B * T_perp_m_cached(q2, bar) - lambda(q2) * (T_perp_m_cached(q2, bar) + T_par_m_cached(q2, bar)) / (mB2 - mK2);
        size_t id = 4 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.a_k_low[id] + cache.b_k_low[id] * q2 / 6.0;
        m_b_local = cache.m_b_PS;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * 8. * cache.m_B * mK2 / (cache.m_B + cache.m_Ks) * cache.ff_calculator.get(BV_FF::T23, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) w += cache.qcdf_calculator.Y(q2);
        F = w * 16. * cache.m_B * mK2 * cache.ff_calculator.get(BV_FF::A12, q2);
        delta_A = delta_A_0(q2, sign, bar);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("----- A_0 -----\n");
    // printf("pref = %.4e + %.4e i\n", real(-N(q2, bar) / (2. * cache.m_Ks * std::sqrt(q2))), imag(-N(q2, bar) / (2. * cache.m_Ks * std::sqrt(q2))));
    // printf("F_V = %.4e + %.4e i\n", real(F), imag(F));
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));
    // printf("delta_A = %.4e + %.4e i\n", real(delta_A), imag(delta_A));

    return (-N(q2, bar) / (2. * cache.m_Ks * std::sqrt(q2)) * (F + 2. * m_b_local * F_T) + delta_A) * had_err_factor;
}

complex_t BKstarllDecay::A_t_low(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoefMapper::cq2_for_lepton_index(static_cast<int>(cfg.gen))] - cache.C[WCoefMapper::cpq2_for_lepton_index(static_cast<int>(cfg.gen))];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }

    complex_t F;
    if (cfg.ff_type == B_FF_Type::SOFT) {
        F = cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PAR, q2) / (cache.m_Ks * cache.qcdf_calculator.Delta_par(q2));
    } else {
        F = cache.ff_calculator.get(BV_FF::A0, q2);
    }

    return N(q2, bar) * std::sqrt(lambda(q2) / q2) * ((2. * C10 + q2 / (cache.m_l * (cache.m_b_mu_b + cache.m_s)) * CQ2)) * F;
}

complex_t BKstarllDecay::A_S_low(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoefMapper::cq1_for_lepton_index(static_cast<int>(cfg.gen))] - cache.C[WCoefMapper::cpq1_for_lepton_index(static_cast<int>(cfg.gen))];
    if (bar) CQ1 = std::conj(CQ1);

    complex_t F;
    if (cfg.ff_type == B_FF_Type::SOFT) {
        F = cache.ff_calculator.E(q2) * cache.ff_calculator.get(BV_FF::XI_PAR, q2) / (cache.m_Ks * cache.qcdf_calculator.Delta_par(q2));
    } else {
        F = cache.ff_calculator.get(BV_FF::A0, q2);
    }

    return -2. * N(q2, bar) * std::sqrt(lambda(q2)) * CQ1 / (cache.m_b_mu_b + cache.m_s) * F;
}

complex_t BKstarllDecay::delta_A_par(double q2, double sign, bool bar) {
    switch(cfg.power_corr_impl) {
        case BKstarllConfig::Power_Corrections_Impl::BFS:
            return delta_A_par_QCDf(q2, sign, bar);
        case BKstarllConfig::Power_Corrections_Impl::BCvDV:
            return delta_A_par_vD(q2, bar);
        case BKstarllConfig::Power_Corrections_Impl::KMPW:
            return delta_A_par_K(q2, bar);
        default:
            return 0.0;
    }
}

complex_t BKstarllDecay::delta_A_0_QCDf(double q2, double sign, bool bar) {
    complex_t guesstimate_err = 1.0;
    complex_t delta_A_PC = 0.0;

    if (!fpeq(std::abs(cache.h_p_fit[0]), 0.0)) {
        complex_t h_0 = cache.h_0_fit[0] + q2 * (cache.h_0_fit[1] + q2 * cache.h_0_fit[2]);
        delta_A_PC = 32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 3) * h_0 / std::sqrt(q2);
    } else {
        size_t id = 4 + size_t (0.5 * (1 + sign));
        guesstimate_err = 1.0 + cache.a_k_low[id] * std::exp(I * cache.phi_k_low[id]) + cache.b_k_low[id] * std::exp(I * cache.theta_k_low[id]) * q2 / 6.0;
    }

    const double mB2 = cache.m_B * cache.m_B;
    const double mB3 = cache.m_B * mB2;
    const double mK2 = cache.m_Ks * cache.m_Ks;
    const double f = lambda(q2) / ((mB2 - mK2) * mB2);

    // Hot path: reuse the lookup/interpolation tables instead of re-triggering the
    // full QCDf calculator inside the binned integrations.
    const complex_t Tperp_m = T_perp_m_cached(q2, bar);
    const complex_t Tpar_m  = T_par_m_cached(q2, bar);

    complex_t delta_A_QCDf =
        -N(q2, bar)
        * (cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI))
        * mB2 / (std::sqrt(q2) * cache.m_Ks)
        * ((2 * (mB2 + 3 * mK2 - q2) * cache.ff_calculator.E(q2) / mB3 - f) * Tperp_m
           - f * Tpar_m);

    return delta_A_QCDf * guesstimate_err + delta_A_PC;
}

complex_t BKstarllDecay::delta_A_0_vD(double q2, bool bar) {
    double F_0  = cache.m_Ks / std::sqrt(q2) * ((std::pow(cache.m_B, 2) - std::pow(cache.m_Ks, 2) - q2) * std::pow(cache.m_B + cache.m_Ks, 2) * cache.ff_calculator.get(BV_FF::A1, q2) - lambda(q2) * cache.ff_calculator.get(BV_FF::A2, q2))/(2. * std::pow(cache.m_Ks, 2) * std::pow(cache.m_B + cache.m_Ks, 2));
    complex_t z_q2 = cache.ff_calculator.z(q2, cache.tp_nf, cache.t0_nf);
    complex_t P_H_0 = cache.alpha_0[0] + z_q2 * cache.alpha_0[1];
    complex_t H_0 = (1. - z_q2 * std::conj(cache.z_Jpsi_nf)) / (z_q2 - cache.z_Jpsi_nf) * (1. - z_q2 * std::conj(cache.z_psi2S_nf)) / (z_q2 - cache.z_psi2S_nf) * P_H_0 * F_0;
    return 32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 2) * (cache.m_B + cache.m_Ks) * H_0 / q2;
}

complex_t BKstarllDecay::delta_A_0_K(double q2, bool bar) {
    double DeltaC9_M2 = (cache.r1_M[1] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[1] * cache.q2_bar / q2) / (1 + cache.r2_M[1] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    double DeltaC9_M3 = (cache.r1_M[2] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[2] * cache.q2_bar / q2) / (1 + cache.r2_M[2] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    return -N(q2, bar) / 2. / cache.m_Ks / std::sqrt(q2) * (((cache.m_B * cache.m_B- cache.m_Ks * cache.m_Ks - q2) * (cache.m_B + cache.m_Ks) * cache.ff_calculator.get(BV_FF::A1, q2) * DeltaC9_M2 - lambda(q2) * cache.ff_calculator.get(BV_FF::A2, q2) * DeltaC9_M3 / (cache.m_B + cache.m_Ks)));
}

complex_t BKstarllDecay::delta_A_0(double q2, double sign, bool bar) {
    switch(cfg.power_corr_impl) {
        case BKstarllConfig::Power_Corrections_Impl::BFS:
            return delta_A_0_QCDf(q2, sign, bar);
        case BKstarllConfig::Power_Corrections_Impl::BCvDV:
            return delta_A_0_vD(q2, bar);
        case BKstarllConfig::Power_Corrections_Impl::KMPW:
            return delta_A_0_K(q2, bar);
        default:
            return 0.0;
    }
}

complex_t BKstarllDecay::C7_eff(double q2, bool bar) {
    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    return (bar ? std::conj(cache.C[WCoef::C7]) : cache.C[WCoef::C7]) + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*BV::f_87(s_hat, cache.L_b));
}

complex_t BKstarllDecay::C9_eff(double q2, bool bar) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+l_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t B = BV::B_Seidel(s_hat, cache.L_b);
    complex_t C = BV::C_Seidel(q2, cache.mu_b);

    return (bar ? std::conj(cache.C[WCoef::C9]) : cache.C[WCoef::C9])
         + BV::h(q2, 0., cache.mu_b) * C_h0
         + BV::h(q2, cache.m_b_PS, cache.mu_b) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * BV::f_89(s_hat))
         + std::pow(cache.m_c_mu_b, 2) / q2 * C_mc;
}

complex_t BKstarllDecay::A_perp_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) + (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) + (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] + cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);

    // printf("C7eff = %.4e + %.4e i\n", std::real(C7_eff(q2, bar)), std::imag(C7_eff(q2, bar)));
    // printf("C9eff = %.4e + %.4e i\n", std::real(C9_eff(q2, bar)), std::imag(C9_eff(q2, bar)));
    // printf("m_B = %.4e\n", cache.m_B);
    // printf("m_Kstar = %.4e\n", cache.m_Ks);
    // printf("N0 = %.4e + %.4e i\n", real(cache.N_0), imag(cache.N_0));
    // printf("N/N0 = %.4e\n", std::real(N(q2, bar) / cache.N_0));
    // printf("N = %.4e + %.4e i\n", std::real(N(q2, bar)), std::imag(N(q2, bar)));
    // printf("f_perp = %.4e\n", cache.ff_calculator.get(BV_FF::F_PERP, q2));

    return N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mb * cache.m_B / q2 * C7) * cache.ff_calculator.get(BV_FF::F_PERP, q2) * (1. + cache.a_k_high[size_t (0.5 * (1 + sign))] * std::exp(I * cache.phi_k_high[size_t (0.5 * (1 + sign))]));
}

complex_t BKstarllDecay::A_par_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mb * cache.m_B / q2 * C7) * cache.ff_calculator.get(BV_FF::F_PAR, q2) * (1. + cache.a_k_high[2 + size_t (0.5 * (1 + sign))] * std::exp(I * cache.phi_k_high[2 + size_t (0.5 * (1 + sign))]));
}

complex_t BKstarllDecay::A_0_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mb * cache.m_B / q2 * C7) * cache.ff_calculator.get(BV_FF::F_0, q2)  * (1. + cache.a_k_high[4 + size_t (0.5 * (1 + sign))] * std::exp(I * cache.phi_k_high[4 + size_t (0.5 * (1 + sign))]));
}

complex_t BKstarllDecay::A_t_high(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoefMapper::cq2_for_lepton_index(static_cast<int>(cfg.gen))] - cache.C[WCoefMapper::cpq2_for_lepton_index(static_cast<int>(cfg.gen))];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }
    return N(q2, bar) * sqrt(lambda(q2) / q2) * (2. * C10 + q2 / cache.m_l * CQ2 / (cache.m_b_mb + cache.m_s)) * cache.ff_calculator.get(BV_FF::A0, q2) * (1. + cache.a_k_high[6] * std::exp(I * cache.phi_k_high[6]));
}

complex_t BKstarllDecay::A_S_high(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoefMapper::cq1_for_lepton_index(static_cast<int>(cfg.gen))] - cache.C[WCoefMapper::cpq1_for_lepton_index(static_cast<int>(cfg.gen))];
    if (bar) CQ1 = std::conj(CQ1);
    return -2. * N(q2, bar) * sqrt(lambda(q2)) * CQ1 / (cache.m_b_mb + cache.m_s) * cache.ff_calculator.get(BV_FF::A0, q2) * (1. + cache.a_k_high[7] * std::exp(I * cache.phi_k_high[7]));
}

complex_t BKstarllDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BKstarllDecay::A_perp(double q2, double sign, bool bar) {
    return interpolate(q2, A_perp_low(q2, sign, bar), A_perp_high(q2, sign, bar));
}

complex_t BKstarllDecay::A_par(double q2, double sign, bool bar) {
    return interpolate(q2, A_par_low(q2, sign, bar), A_par_high(q2, sign, bar));
}

complex_t BKstarllDecay::A_0(double q2, double sign, bool bar) {
    return interpolate(q2, A_0_low(q2, sign, bar), A_0_high(q2, sign, bar));
}

complex_t BKstarllDecay::A_t(double q2, bool bar) {
    return interpolate(q2, A_t_low(q2, bar), A_t_high(q2, bar));
}

complex_t BKstarllDecay::A_S(double q2, bool bar) {
    return interpolate(q2, A_S_low(q2, bar), A_S_high(q2, bar));
}

double BKstarllDecay::J1s(double q2, bool bar) {
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

double BKstarllDecay::J1c(double q2, bool bar) {
    return std::pow(std::abs(A_0(q2, -1, bar)), 2) + std::pow(std::abs(A_0(q2, 1, bar)), 2) 
         + std::pow(2 * cache.m_l, 2) / q2 * (
              std::pow(std::abs(A_t(q2, bar)), 2)
            + 2. * std::real(A_0(q2, -1, bar) * std::conj(A_0(q2, 1, bar)))
           )
         + std::pow(beta_l(q2) * std::abs(A_S(q2, bar)), 2);
}

double BKstarllDecay::J2s(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_par(q2, -1, bar)), 2)
      + std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BKstarllDecay::J2c(double q2, bool bar) {
    return -std::pow(beta_l(q2), 2) * (
        std::pow(std::abs(A_0(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_0(q2, 1, bar)), 2)
    );
}

double BKstarllDecay::J3(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 2. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      - std::pow(std::abs(A_par(q2, -1, bar)), 2)
      - std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BKstarllDecay::J4(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::real(A_0(q2, 1, bar) * std::conj(A_par(q2, 1, bar)))
    );
}

double BKstarllDecay::J5(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
      - cache.m_l / std::sqrt(q2) * std::real((A_par(q2, -1, bar) + A_par(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BKstarllDecay::J6s(double q2, bool bar) {
    return 2. * beta_l(q2) * (
        std::real(A_par(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_par(q2, 1, bar) * std::conj(A_perp(q2, 1, bar))) 
    );
}

double BKstarllDecay::J6c(double q2, bool bar) {
    return 4. * beta_l(q2) * cache.m_l / std::sqrt(q2) * (std::real((A_0(q2, -1, bar) + A_0(q2, 1, bar)) * std::conj(A_S(q2, bar))));
}

double BKstarllDecay::J7(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      - std::imag(A_0(q2, 1, bar) * std::conj(A_par(q2, 1, bar)))
      + cache.m_l / std::sqrt(q2) * std::imag((A_perp(q2, -1, bar) + A_perp(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BKstarllDecay::J8(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      + std::imag(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
    );
}

double BKstarllDecay::J9(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) * (
        std::imag(A_perp(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::imag(A_perp(q2, 1, bar) * std::conj(A_par(q2, 1, bar))) 
    );
}

void BKstarllDecay::compute_binned_J_i() {
    // One pass per bin and per CP-conjugation:
    // - evaluate the 8 transversity amplitudes once per q2 point
    // - derive all J_i from those amplitudes
    // - accumulate every needed binned integral together
    //
    // This removes the previous "15 separate adaptive integrals per bin" pattern,
    // which was recomputing the same amplitudes over and over.
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

    auto eval_integrands = [&] (double q2, bool bar) -> std::array<double, 15> {
        const AmpSet a = eval_amplitudes(q2, bar);
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

        const double j1s =
            (2.0 + a.beta2) / 4.0 * (norm_Aperp_L + norm_Aperp_R + norm_Apar_L + norm_Apar_R)
            + four_ml2_over_q2 * std::real(a.A_perp_L * std::conj(a.A_perp_R)
                                         + a.A_par_L  * std::conj(a.A_par_R));

        const double j1c =
            norm_A0_L + norm_A0_R
            + four_ml2_over_q2 * (norm_At + 2.0 * std::real(a.A_0_L * std::conj(a.A_0_R)))
            + a.beta2 * norm_AS;

        const double j2s = a.beta2 / 4.0 * (norm_Aperp_L + norm_Aperp_R + norm_Apar_L + norm_Apar_R);
        const double j2c = -a.beta2 * (norm_A0_L + norm_A0_R);

        const double j3 = a.beta2 / 2.0 * (norm_Aperp_L + norm_Aperp_R - norm_Apar_L - norm_Apar_R);

        const double j4 = a.beta2 / sqrt2 * (
            std::real(a.A_0_L * std::conj(a.A_par_L))
          + std::real(a.A_0_R * std::conj(a.A_par_R))
        );

        const double j5 = a.beta * sqrt2 * (
            std::real(a.A_0_L * std::conj(a.A_perp_L))
          - std::real(a.A_0_R * std::conj(a.A_perp_R))
          - ml_over_sqrt_q2 * std::real(Apar_sum * std::conj(a.A_S))
        );

        const double j6s = 2.0 * a.beta * (
            std::real(a.A_par_L * std::conj(a.A_perp_L))
          - std::real(a.A_par_R * std::conj(a.A_perp_R))
        );

        const double j6c = 4.0 * a.beta * ml_over_sqrt_q2 * std::real(A0_sum * std::conj(a.A_S));

        const double j7 = a.beta * sqrt2 * (
            std::imag(a.A_0_L * std::conj(a.A_par_L))
          - std::imag(a.A_0_R * std::conj(a.A_par_R))
          + ml_over_sqrt_q2 * std::imag(Aperp_sum * std::conj(a.A_S))
        );

        const double j8 = a.beta2 / sqrt2 * (
            std::imag(a.A_0_L * std::conj(a.A_perp_L))
          + std::imag(a.A_0_R * std::conj(a.A_perp_R))
        );

        const double j9 = a.beta2 * (
            std::imag(a.A_perp_L * std::conj(a.A_par_L))
          + std::imag(a.A_perp_R * std::conj(a.A_par_R))
        );

        return {{
            2.0 * j1s + j1c,
            j2s,
            j2c,
            j3,
            j4,
            j5,
            a.beta * j5,
            j6s,
            a.beta * j6s,
            j6c,
            j7,
            a.beta * j7,
            j8,
            j9,
            j1c
        }};
    };

    auto integrate_bin = [&] (double q2_l, double q2_u, bool bar) -> std::array<double, 15> {
        std::array<double, 15> acc {};
        const double center = 0.5 * (q2_l + q2_u);
        const double half_width = 0.5 * (q2_u - q2_l);

        for (size_t i = 0; i < GL24_X.size(); ++i) {
            const double q2 = center + half_width * GL24_X[i];
            const auto vals = eval_integrands(q2, bar);
            for (size_t k = 0; k < acc.size(); ++k) {
                acc[k] += GL24_W[i] * vals[k];
            }
        }

        for (double& v : acc) {
            v *= half_width;
        }

        return acc;
    };

    const auto& bins = this->bins.value();
    const size_t nbins = bins.size();

    clear_and_reserve(cache.J_i_binned, nbins);
    clear_and_reserve(cache.J_i_bar_binned, nbins);

    auto fill_binned = [&] (std::array<std::vector<double>, 15>& dest, bool bar) {
        for (const auto& [q2_l, q2_u] : bins) {
            const auto integ = integrate_bin(q2_l, q2_u, bar);
            for (size_t k = 0; k < dest.size(); ++k) {
                dest[k].emplace_back(integ[k]);
            }
        }
    };

    fill_binned(cache.J_i_binned, false);
    fill_binned(cache.J_i_bar_binned, true);
}

// std::vector<ObservableValue> BKstarllDecay::dBR_dq2_binned(bool bar, Observables id, bool br) {
//     std::vector<ObservableValue> out;
//     auto J_i = bar ? cache.J_i_bar_binned : cache.J_i_binned;
//     double br_factor = br ? cache.life_B : 1.0;
//     for (size_t i = 0; i < this->bins.value().size(); i++) {
//         double res = 0.75 * (cache.J_i_binned[0][i] + cache.J_i_bar_binned[0][i] - (2 * (cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i]) + cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i]) / 3.) * br_factor; 
//         out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
//     }   
//     return out;
// }

std::vector<ObservableValue> BKstarllDecay::dBR_dq2_binned(bool bar, Observables id, bool br) {
    (void)bar; // current public observables are CP-averaged: B + Bbar, then /2 below.

    std::vector<ObservableValue> out;
    const auto& bins = this->bins.value();
    const double br_factor = br ? cache.life_B : 1.0;

    for (size_t i = 0; i < bins.size(); i++) {
        const auto& [q2_l, q2_u] = bins[i];
        const double width = q2_u - q2_l;

        if (!std::isfinite(width) || width <= 0.0) {
            LOG_WARN("BKstarll dBR/dq2: invalid bin width", q2_l, q2_u);
            out.emplace_back(ObservableMapper::to_id(id),
                             std::numeric_limits<double>::quiet_NaN(),
                             bins[i]);
            continue;
        }

        const double gamma_integrated_sum =
            0.75 * (
                cache.J_i_binned[0][i] + cache.J_i_bar_binned[0][i]
                - (
                    2.0 * (cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i])
                    + cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i]
                  ) / 3.0
            );

        const double res = 0.5 * gamma_integrated_sum * br_factor / width;
        out.emplace_back(ObservableMapper::to_id(id), res, bins[i]);
    }

    return out;
}

double BKstarllDecay::dG_dq2_avg_bin(size_t bin) {
    return 0.75 * (cache.J_i_binned[0][bin] + cache.J_i_bar_binned[0][bin] - (2 * (cache.J_i_binned[1][bin] + cache.J_i_bar_binned[1][bin]) + cache.J_i_binned[2][bin] + cache.J_i_bar_binned[2][bin]) / 3.);
}

std::vector<ObservableValue> BKstarllDecay::A_FB_binned(Observables id, bool cpv) {
    std::vector<ObservableValue> out;
    double sign = cpv ? -1 : 1;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J6 = 2 * cache.J_i_binned[7][i] + sign * cache.J_i_binned[9][i];
        double J6bar = 2 * cache.J_i_bar_binned[7][i] + sign * cache.J_i_bar_binned[9][i];
        double res = -0.375 * (J6 + J6bar) / dG_dq2_avg_bin(i); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

ObservableValue BKstarllDecay::q0(Observables id) {
    auto f = [this] (double q2) {
        return 2 * (J6s(q2, true) + J6s(q2, false)) + J6c(q2, true) + J6c(q2, false); 
    };

    double q2a, q2b;
    bool found_bracket = find_bracket(f, 1.01 * cache.q2_min, 0.99 * cache.q2_max, q2a, q2b);
    if (!found_bracket) {
        LOG_WARN("Forwards-Backwards asymmetry in B > K*ll doesn't cross 0.");
        return ObservableValue(ObservableMapper::to_id(id), NAN);
    } 

    return ObservableValue(ObservableMapper::to_id(id), brent_root(f, q2a, q2b));
}

std::vector<ObservableValue> BKstarllDecay::A_CP_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double dG = cache.J_i_binned[0][i] - (2 * cache.J_i_binned[1][i] + cache.J_i_binned[2][i]) / 3.;
        double dGbar = cache.J_i_bar_binned[0][i] - (2 * cache.J_i_bar_binned[1][i] + cache.J_i_bar_binned[2][i]) / 3.;
        double res = (dG - dGbar) / (dG + dGbar); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::F_L_binned(Observables id) {
    std::vector<ObservableValue> out;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        const double num =
            0.75 * (cache.J_i_binned[14][i] + cache.J_i_bar_binned[14][i])
            - 0.25 * (cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i]);

        const double den = dG_dq2_avg_bin(i);

        const double res = safe_div(num, den, "BKstarll F_L");

        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }

    return out;
}

std::vector<ObservableValue> BKstarllDecay::F_T_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = 4.0 * (cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i]) / dG_dq2_avg_bin(i); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_1_binned(Observables id) {
    auto num_f = [this] (double q2) {
        double AperpApar = std::real(A_par(q2, 1, false) * std::conj(A_perp(q2, 1, false)) + A_par(q2, -1, false) * std::conj(A_perp(q2, -1, false)));
        double AperpApar_bar = std::real(A_par(q2, 1, true) * std::conj(A_perp(q2, 1, true)) + A_par(q2, -1, true) * std::conj(A_perp(q2, -1, true)));
        return std::pow(beta_l(q2), 2) * std::real(AperpApar + AperpApar_bar);
    };

    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scpa = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double num = integrate(num_f, this->bins.value()[i].first, this->bins.value()[i].second, 1e-2);
        double res = -0.5 * num / J2scpa; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_2_binned(Observables id) {
    std::vector<ObservableValue> out;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        const double num = 0.5 * (cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i]);
        const double den = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];

        const double res = safe_div(num, den, "BKstarll A_T_2");

        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }

    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_3_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J4cp = cache.J_i_binned[4][i] + cache.J_i_bar_binned[4][i];
        double J7cp = cache.J_i_binned[11][i] + cache.J_i_bar_binned[11][i];
        double res = std::sqrt((4 * J4cp * J4cp + J7cp * J7cp) / std::abs(-2 * J2ccp * (2 * J2scp + J3cp)));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_4_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J4cp = cache.J_i_binned[4][i] + cache.J_i_bar_binned[4][i];
        double J5cp = cache.J_i_binned[6][i] + cache.J_i_bar_binned[6][i];
        double J7cp = cache.J_i_binned[11][i] + cache.J_i_bar_binned[11][i];
        double J8cp = cache.J_i_binned[12][i] + cache.J_i_bar_binned[12][i];
        double res = std::sqrt((J5cp * J5cp + 4 * J8cp * J8cp) / (J7cp * J7cp + 4 * J4cp * J4cp));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_5_binned(Observables id) {
    auto num_f = [this] (double q2) {
        complex_t AperpApar = A_par(q2, 1, false) * std::conj(A_perp(q2, -1, false)) + A_perp(q2, 1, false) * std::conj(A_par(q2, -1, false));
        complex_t AperpApar_bar = A_par(q2, 1, true) * std::conj(A_perp(q2, -1, true)) + A_perp(q2, 1, true) * std::conj(A_par(q2, -1, true));
        return std::abs(AperpApar + AperpApar_bar);
    };

    auto den_f = [this] (double q2) {
        double Aperp2Apar2 = (std::pow(std::abs(A_perp(q2, -1, false)), 2) + std::pow(std::abs(A_perp(q2, 1, false)), 2) + std::pow(std::abs(A_par(q2, -1, false)), 2) + std::pow(std::abs(A_par(q2, 1, false)), 2));
        double Aperp2Apar2bar = (std::pow(std::abs(A_perp(q2, -1, true)), 2) + std::pow(std::abs(A_perp(q2, 1, true)), 2) + std::pow(std::abs(A_par(q2, -1, true)), 2) + std::pow(std::abs(A_par(q2, 1, true)), 2));
        return Aperp2Apar2 + Aperp2Apar2bar;
    };

    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double num = integrate(num_f, this->bins.value()[i].first, this->bins.value()[i].second, 1e-2);
        double den = integrate(den_f, this->bins.value()[i].first, this->bins.value()[i].second, 1e-2);
        double res = num / den; 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_Re_binned(Observables id) {
    std::vector<ObservableValue> out;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        const double num = 0.25 * (cache.J_i_binned[8][i] + cache.J_i_bar_binned[8][i]);
        const double den = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];

        const double res = safe_div(num, den, "BKstarll A_T_Re");

        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }

    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_Re_CPV_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J6scpa = cache.J_i_binned[8][i] - cache.J_i_bar_binned[8][i];
        double res = 0.25 * J6scpa / J2scp;
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_Im_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = (cache.J_i_binned[13][i] + cache.J_i_bar_binned[13][i]) / dG_dq2_avg_bin(i); 
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::alpha_K_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double res = -0.5 * (2 * J2scp + J2ccp) / J2scp;
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::H_T_1_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J4cp = cache.J_i_binned[4][i] + cache.J_i_bar_binned[4][i];
        double res = RT2 * J4cp / std::sqrt(std::abs(J2ccp * (J2scp - J3cp)));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::H_T_2_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J5cp = cache.J_i_binned[6][i] + cache.J_i_bar_binned[6][i];
        double res = J5cp / std::sqrt(std::abs(2 * J2ccp * (2 * J2scp + J3cp)));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::H_T_3_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J6cp = (2 * (cache.J_i_binned[7][i] + cache.J_i_bar_binned[7][i]) + cache.J_i_binned[9][i] + cache.J_i_bar_binned[9][i]);
        double res = 0.5 * J6cp / std::sqrt(std::abs(4 * J2scp * J2scp - J3cp * J3cp));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_2_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J6scp = cache.J_i_binned[7][i] + cache.J_i_bar_binned[7][i];
        double res = 0.125 * J6scp / J2scp;
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_3_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J9cp = cache.J_i_binned[13][i] + cache.J_i_bar_binned[13][i];
        double res = -0.25 * J9cp / J2scp;
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_6_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J7cp = cache.J_i_binned[11][i] + cache.J_i_bar_binned[11][i];
        double res = -J7cp / std::sqrt(std::abs(2 * J2ccp * (2 * J2scp - J3cp)));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_8_binned(Observables id) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J8cp = cache.J_i_binned[12][i] + cache.J_i_bar_binned[12][i];
        double res = -RT2 * J8cp / std::sqrt(std::abs(J2ccp * (2 * J2scp - J3cp)));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::Pp_i_binned(size_t i, bool cpv, Observables id) {
    if (!(i == 4 || i == 5 || i == 6 || i == 8)) LOG_ERROR("Value Error", "P'_i(B > K*ll) is not defined for i =", i);

    std::map<size_t, double> factors = {{4, 1.0}, {5, 0.5}, {6, -0.5}, {8, -1.0}};
    std::map<size_t, size_t> J_idx = {{4, 4}, {5, 5}, {6, 10}, {8, 12}};
    double sign = cpv ? -1 : 1;

    std::vector<ObservableValue> out;
    for (size_t j = 0; j < this->bins.value().size(); j++) {
        double J2scp = cache.J_i_binned[1][j] + cache.J_i_bar_binned[1][j];
        double J2ccp = cache.J_i_binned[2][j] + cache.J_i_bar_binned[2][j];
        double Jicp = cache.J_i_binned[J_idx[i]][j] + sign * cache.J_i_bar_binned[J_idx[i]][j];
        double res = factors[i] * Jicp / std::sqrt(std::abs(J2ccp * J2scp));
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[j]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::S_i_binned(size_t i, bool cpv, Observables id) {
    if (i > 10) LOG_ERROR("Value Error", "S_i(B > K*ll) is not defined for i =", i);

    // 0: S1c, 1: S2s, 2: S2c
    std::map<size_t, size_t> J_idx = {{0, 14}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 7}, {7, 10}, {8, 12}, {9, 13}, {10, 9}};
    double sign = cpv ? -1 : 1;
    double global_sign = i == 10 || i == 2 && cpv ? -1 : 1;

    std::vector<ObservableValue> out;
    for (size_t j = 0; j < this->bins.value().size(); j++) {
        double res = global_sign * (cache.J_i_binned[J_idx[i]][j] + sign * cache.J_i_bar_binned[J_idx[i]][j]) / dG_dq2_avg_bin(j); 
        out.emplace_back(ObservableMapper::to_id(id), i == 6 ? -res : res, this->bins.value()[j]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_i_CPV_binned(size_t i, Observables id) {
   if (i < 1 || i > 3) LOG_ERROR("Value Error", "P_i_CPV(B > K*ll) is not defined for i =", i);

    std::map<size_t, double> factors = {{1, 0.5}, {2, 0.125}, {3, -0.25}};
    std::map<size_t, size_t> J_idx = {{1, 3}, {2, 7}, {3, 13}};

    std::vector<ObservableValue> out;
    for (size_t j = 0; j < this->bins.value().size(); j++) {
        double J2scp = cache.J_i_binned[1][j] + cache.J_i_bar_binned[1][j];
        double Jicpv = cache.J_i_binned[J_idx[i]][j] - cache.J_i_bar_binned[J_idx[i]][j];
        double res = factors[i] * Jicpv / J2scp;
        out.emplace_back(ObservableMapper::to_id(id), res, this->bins.value()[j]);
    }   
    return out;
}

void BKstarllDecay::test_ff() {
    std::ofstream fs;
    fs.open("B_Ksll_FF.csv");
    fs << "q2,A0,A1,A12,V,T1,T2,T23\n";

    auto write_line = [&] (double q2) {
        fs << q2 << "," 
        << cache.ff_calculator.get(BV_FF::A0, q2) << ","
        << cache.ff_calculator.get(BV_FF::A1, q2) << ","
        << cache.ff_calculator.get(BV_FF::A12, q2) << ","
        << cache.ff_calculator.get(BV_FF::V, q2) << ","
        << cache.ff_calculator.get(BV_FF::T1, q2) << ","
        << cache.ff_calculator.get(BV_FF::T2, q2) << ","
        << cache.ff_calculator.get(BV_FF::T23, q2)
        << "\n";
    };

    double q2_min = 4 * std::pow(0.1057, 2);
    double q2_max = 19.2542;
    size_t n = 200;
    double dq2 = (q2_max - q2_min) / n;
    double q2 = q2_min;
    for (size_t i = 0; i <= n; i++) {
        write_line(q2);
        q2 += dq2;
    }
}

void BKstarllDecay::test_T() {
    std::ofstream fs;
    fs.open("B_Ksll_T.csv");
    fs << "q2,T_perp_p_re,T_perp_p_im,T_perp_m_re,T_perp_m_im,T_par_p_re,T_par_p_im,T_par_m_re,T_par_m_im\n";

    auto write_line = [&] (double q2) {
        fs << q2 
        << "," << std::real(T_perp_p_cached(q2, false)) << "," << std::imag(T_perp_p_cached(q2, false))
        << "," << std::real(T_perp_m_cached(q2, false)) << "," << std::imag(T_perp_m_cached(q2, false))
        << "," << std::real(T_par_m_cached(q2, false)) << "," << std::imag(T_par_m_cached(q2, false))
        << "\n";
    };

    size_t n = 200;
    double dq2 = (cache.q2_high - cache.q2_min) / n;
    double q2 = cache.q2_min;
    for (size_t i = 0; i <= n; i++) {
        write_line(q2);
        q2 += dq2;
    }
}

void BKstarllDecay::test_J() {
    std::ofstream fs;
    fs.open("B_Ksll_J.csv");
    fs << "q2,J1s,J1c,J2s,J2c,J3,J4,J5,J6s,J6c,J7,J8,J9,J1sbar,J1cbar,J2sbar,J2cbar,J3bar,J4bar,J5bar,J6sbar,J6cbar,J7bar,J8bar,J9bar\n";

    auto write_line = [&] (double q2) {
        fs << q2 
        << "," << J1s(q2, false)
        << "," << J1c(q2, false)
        << "," << J2s(q2, false)
        << "," << J2c(q2, false)
        << "," << J3(q2, false)
        << "," << J4(q2, false)
        << "," << J5(q2, false)
        << "," << J6s(q2, false)
        << "," << J6c(q2, false)
        << "," << J7(q2, false)
        << "," << J8(q2, false)
        << "," << J9(q2, false)
        << "," << J1s(q2, true)
        << "," << J1c(q2, true)
        << "," << J2s(q2, true)
        << "," << J2c(q2, true)
        << "," << J3(q2, true)
        << "," << J4(q2, true)
        << "," << J5(q2, true)
        << "," << J6s(q2, true)
        << "," << J6c(q2, true)
        << "," << J7(q2, true)
        << "," << J8(q2, true)
        << "," << J9(q2, true)
        << "\n";
    };

    size_t n = 200;
    double dq2 = (cache.q2_max - cache.q2_min) / n;
    double q2 = cache.q2_min;

    for (size_t i = 0; i <= n; i++) {
        write_line(q2);
        q2 += dq2;
    }
}

std::vector<ObservableValue> BKstarllDecay::Rm1_BKstar(
    Observables id,
    BKstarllConfig::B_Charge charge
) {
    std::vector<ObservableValue> out;
    std::vector<double> Gamma_mu;
    std::vector<double> Gamma_e;

    set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, charge);

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        Gamma_mu.emplace_back(dG_dq2_avg_bin(i));
    }

    set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, charge);

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        Gamma_e.emplace_back(dG_dq2_avg_bin(i));
    }

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        out.emplace_back(
            ObservableMapper::to_id(id),
            Gamma_mu[i] / Gamma_e[i] - 1.0,
            this->bins.value()[i]
        );
    }

    return out;
}


std::vector<ObservableValue> BKstarllDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::DBR_DQ2_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return dBR_dq2_binned(false, obs, true);   
    case Observables::DGAMMA_DQ2_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return dBR_dq2_binned(false, obs, false);   
    case Observables::A_FB_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return A_FB_binned(obs, false);
    case Observables::A_FB_CPV_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return A_FB_binned(obs, true);
    case Observables::Q0_A_FB_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return {q0(obs)};
    case Observables::A_CP_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return A_CP_binned(obs);
    case Observables::F_L_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return F_L_binned(obs);
    case Observables::F_T_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return F_T_binned(obs);
    case Observables::A_T_1_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_1_binned(obs);
    case Observables::A_T_2_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_2_binned(obs);
    case Observables::A_T_3_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_3_binned(obs);
    case Observables::A_T_4_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return A_T_4_binned(obs);
    case Observables::A_T_5_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return A_T_5_binned(obs);
    case Observables::A_T_RE_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return A_T_Re_binned(obs);
    case Observables::A_T_RE_CPV_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return A_T_Re_CPV_binned(obs);
    case Observables::A_IM_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return A_Im_binned(obs);
    case Observables::ALPHA_K_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return alpha_K_binned(obs);
    case Observables::H_T_1_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return H_T_1_binned(obs);
    case Observables::H_T_2_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return H_T_2_binned(obs);
    case Observables::H_T_3_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return H_T_3_binned(obs);
    case Observables::P_1_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return A_T_2_binned(obs);
    case Observables::P_2_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return P_2_binned(obs);
    case Observables::P_3_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return P_3_binned(obs);
    case Observables::P_4_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return H_T_1_binned(obs);
    case Observables::P_5_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return H_T_2_binned(obs);
    case Observables::P_6_B__KSTAR_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);   
        return P_6_binned(obs);
    case Observables::P_8_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return P_8_binned(obs);
    case Observables::P_PRIME_4_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(4, false, obs);
    case Observables::P_PRIME_5_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(5, false, obs);
    case Observables::P_PRIME_6_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(6, false, obs);
    case Observables::P_PRIME_8_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(8, false, obs);
    case Observables::S_1C_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(0, false, obs);
    case Observables::S_2S_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(1, false, obs);
    case Observables::S_2C_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(2, false, obs);
    case Observables::S_3_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(3, false, obs);
    case Observables::S_4_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(4, false, obs);
    case Observables::S_5_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(5, false, obs);
    case Observables::S_6C_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return S_i_binned(10, false, obs);
    case Observables::S_7_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(7, false, obs);
    case Observables::S_8_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(8, false, obs);
    case Observables::S_9_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(9, false, obs);
    case Observables::A_1C_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(0, true, obs);
    case Observables::A_2S_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(1, true, obs);
    case Observables::A_FL_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(2, true, obs);
    case Observables::A_3_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(3, true, obs);
    case Observables::A_4_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(4, true, obs);
    case Observables::A_5_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(5, true, obs);
    case Observables::A_6S_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(6, true, obs);
    case Observables::A_6C_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(10, true, obs);
    case Observables::A_7_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(7, true, obs);
    case Observables::A_8_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(8, true, obs);
    case Observables::A_9_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(9, true, obs);
    case Observables::P_1_CPV_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return P_i_CPV_binned(1, obs);
    case Observables::P_2_CPV_B__KSTAR_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);
        return P_i_CPV_binned(2, obs);
    case Observables::P_3_CPV_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return P_i_CPV_binned(3, obs);
    case Observables::P_PRIME_4_CPV_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(4, true, obs);
    case Observables::P_PRIME_5_CPV_B__KSTAR_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(5, true, obs);
    case Observables::P_PRIME_6_CPV_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(6, true, obs);
    case Observables::P_PRIME_8_CPV_B__KSTAR_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(8, true, obs);
    case Observables::DBR_DQ2_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return dBR_dq2_binned(false, obs);  
    case Observables::DGAMMA_DQ2_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return dBR_dq2_binned(false, obs, false);    
    case Observables::A_FB_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return A_FB_binned(obs, false);
    case Observables::A_FB_CPV_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return A_FB_binned(obs, true);
    case Observables::Q0_A_FB_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return {q0(obs)};
    case Observables::A_CP_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return A_CP_binned(obs);
    case Observables::F_L_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return F_L_binned(obs);
    case Observables::F_T_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return F_T_binned(obs);
    case Observables::A_T_1_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_1_binned(obs);
    case Observables::A_T_2_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_2_binned(obs);
    case Observables::A_T_3_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_3_binned(obs);
    case Observables::A_T_4_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return A_T_4_binned(obs);
    case Observables::A_T_5_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return A_T_5_binned(obs);
    case Observables::A_T_RE_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return A_T_Re_binned(obs);
    case Observables::A_T_RE_CPV_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return A_T_Re_CPV_binned(obs);
    case Observables::A_IM_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return A_Im_binned(obs);
    case Observables::ALPHA_K_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return alpha_K_binned(obs);
    case Observables::H_T_1_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return H_T_1_binned(obs);
    case Observables::H_T_2_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return H_T_2_binned(obs);
    case Observables::H_T_3_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return H_T_3_binned(obs);
    case Observables::P_1_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return A_T_2_binned(obs);
    case Observables::P_2_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return P_2_binned(obs);
    case Observables::P_3_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return P_3_binned(obs);
    case Observables::P_4_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return H_T_1_binned(obs);
    case Observables::P_5_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return H_T_2_binned(obs);
    case Observables::P_6_B__KSTAR_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);   
        return P_6_binned(obs);
    case Observables::P_8_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return P_8_binned(obs);
    case Observables::P_PRIME_4_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(4, false, obs);
    case Observables::P_PRIME_5_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(5, false, obs);
    case Observables::P_PRIME_6_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(6, false, obs);
    case Observables::P_PRIME_8_B__KSTAR_MU_MU:    
        // LOG_INFO("Computing P'8(B > K* mu mu)"); 
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(8, false, obs);
    case Observables::S_1C_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(0, false, obs);
    case Observables::S_2S_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(1, false, obs);
    case Observables::S_2C_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(2, false, obs);
    case Observables::S_3_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(3, false, obs);
    case Observables::S_4_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(4, false, obs);
    case Observables::S_5_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(5, false, obs);
    case Observables::S_6C_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return S_i_binned(10, false, obs);
    case Observables::S_7_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(7, false, obs);
    case Observables::S_8_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(8, false, obs);
    case Observables::S_9_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(9, false, obs);
    case Observables::A_1C_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(0, true, obs);
    case Observables::A_2S_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(1, true, obs);
    case Observables::A_FL_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(2, true, obs);
    case Observables::A_3_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(3, true, obs);
    case Observables::A_4_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(4, true, obs);
    case Observables::A_5_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(5, true, obs);
    case Observables::A_6S_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(6, true, obs);
    case Observables::A_6C_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(10, true, obs);
    case Observables::A_7_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(7, true, obs);
    case Observables::A_8_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(8, true, obs);
    case Observables::A_9_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(9, true, obs);
    case Observables::P_1_CPV_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return P_i_CPV_binned(1, obs);
    case Observables::P_2_CPV_B__KSTAR_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);
        return P_i_CPV_binned(2, obs);
    case Observables::P_3_CPV_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return P_i_CPV_binned(3, obs);
    case Observables::P_PRIME_4_CPV_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(4, true, obs);
    case Observables::P_PRIME_5_CPV_B__KSTAR_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(5, true, obs);
    case Observables::P_PRIME_6_CPV_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(6, true, obs);
    case Observables::P_PRIME_8_CPV_B__KSTAR_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(8, true, obs);
    case Observables::DBR_DQ2_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return dBR_dq2_binned(false, obs);   
    case Observables::DGAMMA_DQ2_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return dBR_dq2_binned(false, obs, false);   
    case Observables::A_FB_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return A_FB_binned(obs, false);
    case Observables::A_FB_CPV_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return A_FB_binned(obs, true);
    case Observables::Q0_A_FB_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return {q0(obs)};
    case Observables::A_CP_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return A_CP_binned(obs);
    case Observables::F_L_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return F_L_binned(obs);
    case Observables::F_T_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return F_T_binned(obs);
    case Observables::A_T_1_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_1_binned(obs);
    case Observables::A_T_2_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_2_binned(obs);
    case Observables::A_T_3_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return A_T_3_binned(obs);
    case Observables::A_T_4_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return A_T_4_binned(obs);
    case Observables::A_T_5_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return A_T_5_binned(obs);
    case Observables::A_T_RE_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return A_T_Re_binned(obs);
    case Observables::A_T_RE_CPV_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return A_T_Re_CPV_binned(obs);
    case Observables::A_IM_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return A_Im_binned(obs);
    case Observables::ALPHA_K_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return alpha_K_binned(obs);
    case Observables::H_T_1_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return H_T_1_binned(obs);
    case Observables::H_T_2_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return H_T_2_binned(obs);
    case Observables::H_T_3_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return H_T_3_binned(obs);
    case Observables::P_1_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return A_T_2_binned(obs);
    case Observables::P_2_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return P_2_binned(obs);
    case Observables::P_3_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return P_3_binned(obs);
    case Observables::P_4_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return H_T_1_binned(obs);
    case Observables::P_5_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return H_T_2_binned(obs);
    case Observables::P_6_B__KSTAR_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);   
        return P_6_binned(obs);
    case Observables::P_8_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return P_8_binned(obs);
    case Observables::P_PRIME_4_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(4, false, obs);
    case Observables::P_PRIME_5_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(5, false, obs);
    case Observables::P_PRIME_6_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(6, false, obs);
    case Observables::P_PRIME_8_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(8, false, obs);
    case Observables::S_1C_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(0, false, obs);
    case Observables::S_2S_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(1, false, obs);
    case Observables::S_2C_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(2, false, obs);
    case Observables::S_3_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(3, false, obs);
    case Observables::S_4_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(4, false, obs);
    case Observables::S_5_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(5, false, obs);
    case Observables::S_6C_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return S_i_binned(10, false, obs);
    case Observables::S_7_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(7, false, obs);
    case Observables::S_8_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(8, false, obs);
    case Observables::S_9_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(9, false, obs);
    case Observables::A_1C_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(0, true, obs);
    case Observables::A_2S_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(1, true, obs);
    case Observables::A_FL_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(2, true, obs);
    case Observables::A_3_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return S_i_binned(3, true, obs);
    case Observables::A_4_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(4, true, obs);
    case Observables::A_5_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(5, true, obs);
    case Observables::A_6S_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(6, true, obs);
    case Observables::A_6C_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(10, true, obs);
    case Observables::A_7_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(7, true, obs);
    case Observables::A_8_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(8, true, obs);
    case Observables::A_9_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return S_i_binned(9, true, obs);
    case Observables::P_1_CPV_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return P_i_CPV_binned(1, obs);
    case Observables::P_2_CPV_B__KSTAR_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);
        return P_i_CPV_binned(2, obs);
    case Observables::P_3_CPV_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return P_i_CPV_binned(3, obs);
    case Observables::P_PRIME_4_CPV_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(4, true, obs);
    case Observables::P_PRIME_5_CPV_B__KSTAR_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS);  
        return Pp_i_binned(5, true, obs);
    case Observables::P_PRIME_6_CPV_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(6, true, obs);
    case Observables::P_PRIME_8_CPV_B__KSTAR_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_PLUS); 
        return Pp_i_binned(8, true, obs);
    case Observables::DBR_DQ2_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return dBR_dq2_binned(false, obs);   
    case Observables::DGAMMA_DQ2_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return dBR_dq2_binned(false, obs, false);   
    case Observables::A_FB_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return A_FB_binned(obs, false);  
    case Observables::A_FB_CPV_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return A_FB_binned(obs, true);
    case Observables::Q0_A_FB_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return {q0(obs)};
    case Observables::A_CP_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return A_CP_binned(obs);
    case Observables::F_L_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return F_L_binned(obs);
    case Observables::F_T_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return F_T_binned(obs);
    case Observables::A_T_1_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return A_T_1_binned(obs);
    case Observables::A_T_2_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return A_T_2_binned(obs);
    case Observables::A_T_3_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return A_T_3_binned(obs);
    case Observables::A_T_4_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return A_T_4_binned(obs);
    case Observables::A_T_5_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return A_T_5_binned(obs);
    case Observables::A_T_RE_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return A_T_Re_binned(obs);
    case Observables::A_T_RE_CPV_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return A_T_Re_CPV_binned(obs);
    case Observables::A_IM_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return A_Im_binned(obs);
    case Observables::ALPHA_K_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return alpha_K_binned(obs);
    case Observables::H_T_1_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return H_T_1_binned(obs);
    case Observables::H_T_2_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return H_T_2_binned(obs);
    case Observables::H_T_3_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return H_T_3_binned(obs);
    case Observables::P_1_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return A_T_2_binned(obs);
    case Observables::P_2_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return P_2_binned(obs);
    case Observables::P_3_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return P_3_binned(obs);
    case Observables::P_4_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return H_T_1_binned(obs);
    case Observables::P_5_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return H_T_2_binned(obs);
    case Observables::P_6_B0__KSTAR0_E_E:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);   
        return P_6_binned(obs);
    case Observables::P_8_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return P_8_binned(obs);
    case Observables::P_PRIME_4_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(4, false, obs);
    case Observables::P_PRIME_5_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(5, false, obs);
    case Observables::P_PRIME_6_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(6, false, obs);
    case Observables::P_PRIME_8_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(8, false, obs);
    case Observables::S_1C_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(0, false, obs);
    case Observables::S_2S_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(1, false, obs);
    case Observables::S_2C_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(2, false, obs);
    case Observables::S_3_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(3, false, obs);
    case Observables::S_4_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(4, false, obs);
    case Observables::S_5_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(5, false, obs);
    case Observables::S_6C_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return S_i_binned(10, false, obs);
    case Observables::S_7_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(7, false, obs);
    case Observables::S_8_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(8, false, obs);
    case Observables::S_9_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(9, false, obs);
    case Observables::A_1C_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(0, true, obs);
    case Observables::A_2S_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(1, true, obs);
    case Observables::A_FL_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(2, true, obs);
    case Observables::A_3_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(3, true, obs);
    case Observables::A_4_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(4, true, obs);
    case Observables::A_5_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(5, true, obs);
    case Observables::A_6S_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(6, true, obs);
    case Observables::A_6C_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(10, true, obs);
    case Observables::A_7_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(7, true, obs);
    case Observables::A_8_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(8, true, obs);
    case Observables::A_9_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(9, true, obs);
    case Observables::P_1_CPV_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return P_i_CPV_binned(1, obs);
    case Observables::P_2_CPV_B0__KSTAR0_E_E:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);
        return P_i_CPV_binned(2, obs);
    case Observables::P_3_CPV_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return P_i_CPV_binned(3, obs);
    case Observables::P_PRIME_4_CPV_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(4, true, obs);
    case Observables::P_PRIME_5_CPV_B0__KSTAR0_E_E:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(5, true, obs);
    case Observables::P_PRIME_6_CPV_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(6, true, obs);
    case Observables::P_PRIME_8_CPV_B0__KSTAR0_E_E:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::E, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(8, true, obs);
    case Observables::DBR_DQ2_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return dBR_dq2_binned(false, obs);   
    case Observables::DGAMMA_DQ2_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return dBR_dq2_binned(false, obs, false);   
    case Observables::A_FB_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return A_FB_binned(obs, false); 
    case Observables::A_FB_CPV_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return A_FB_binned(obs, true);
    case Observables::Q0_A_FB_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return {q0(obs)};
    case Observables::A_CP_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return A_CP_binned(obs);
    case Observables::F_L_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return F_L_binned(obs);
    case Observables::F_T_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return F_T_binned(obs);
    case Observables::A_T_1_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return A_T_1_binned(obs);
    case Observables::A_T_2_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return A_T_2_binned(obs);
    case Observables::A_T_3_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return A_T_3_binned(obs);
    case Observables::A_T_4_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return A_T_4_binned(obs);
    case Observables::A_T_5_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return A_T_5_binned(obs);
    case Observables::A_T_RE_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return A_T_Re_binned(obs);
    case Observables::A_T_RE_CPV_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return A_T_Re_CPV_binned(obs);
    case Observables::A_IM_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return A_Im_binned(obs);
    case Observables::ALPHA_K_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return alpha_K_binned(obs);
    case Observables::H_T_1_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return H_T_1_binned(obs);
    case Observables::H_T_2_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return H_T_2_binned(obs);
    case Observables::H_T_3_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return H_T_3_binned(obs);
    case Observables::P_1_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return A_T_2_binned(obs);
    case Observables::P_2_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return P_2_binned(obs);
    case Observables::P_3_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return P_3_binned(obs);
    case Observables::P_4_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return H_T_1_binned(obs);
    case Observables::P_5_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return H_T_2_binned(obs);
    case Observables::P_6_B0__KSTAR0_MU_MU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);   
        return P_6_binned(obs);
    case Observables::P_8_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return P_8_binned(obs);
    case Observables::P_PRIME_4_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(4, false, obs);
    case Observables::P_PRIME_5_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(5, false, obs);
    case Observables::P_PRIME_6_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(6, false, obs);
    case Observables::P_PRIME_8_B0__KSTAR0_MU_MU:  
        // LOG_INFO("Computing P'8(B0 > K*0 mu mu)");
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(8, false, obs);
    case Observables::S_1C_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(0, false, obs);
    case Observables::S_2S_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(1, false, obs);
    case Observables::S_2C_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(2, false, obs);
    case Observables::S_3_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(3, false, obs);
    case Observables::S_4_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(4, false, obs);
    case Observables::S_5_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(5, false, obs);
    case Observables::S_6C_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return S_i_binned(10, false, obs);
    case Observables::S_7_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(7, false, obs);
    case Observables::S_8_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(8, false, obs);
    case Observables::S_9_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(9, false, obs);
    case Observables::A_1C_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(0, true, obs);
    case Observables::A_2S_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(1, true, obs);
    case Observables::A_FL_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(2, true, obs);
    case Observables::A_3_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(3, true, obs);
    case Observables::A_4_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(4, true, obs);
    case Observables::A_5_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(5, true, obs);
    case Observables::A_6S_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(6, true, obs);
    case Observables::A_6C_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(10, true, obs);
    case Observables::A_7_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(7, true, obs);
    case Observables::A_8_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(8, true, obs);
    case Observables::A_9_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(9, true, obs);
    case Observables::P_1_CPV_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return P_i_CPV_binned(1, obs);
    case Observables::P_2_CPV_B0__KSTAR0_MU_MU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);
        return P_i_CPV_binned(2, obs);
    case Observables::P_3_CPV_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return P_i_CPV_binned(3, obs);
    case Observables::P_PRIME_4_CPV_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(4, true, obs);
    case Observables::P_PRIME_5_CPV_B0__KSTAR0_MU_MU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(5, true, obs);
    case Observables::P_PRIME_6_CPV_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(6, true, obs);
    case Observables::P_PRIME_8_CPV_B0__KSTAR0_MU_MU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::MU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(8, true, obs);
    case Observables::DBR_DQ2_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return dBR_dq2_binned(false, obs);  
    case Observables::DGAMMA_DQ2_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return dBR_dq2_binned(false, obs, false);    
    case Observables::A_FB_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return A_FB_binned(obs, false);
    case Observables::A_FB_CPV_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return A_FB_binned(obs, true);
    case Observables::Q0_A_FB_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return {q0(obs)};
    case Observables::A_CP_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return A_CP_binned(obs);
    case Observables::F_L_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return F_L_binned(obs);
    case Observables::F_T_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return F_T_binned(obs);
    case Observables::A_T_1_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);   
        return A_T_1_binned(obs);
    case Observables::A_T_2_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);   
        return A_T_2_binned(obs);
    case Observables::A_T_3_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);   
        return A_T_3_binned(obs);
    case Observables::A_T_4_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return A_T_4_binned(obs);
    case Observables::A_T_5_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return A_T_5_binned(obs);
    case Observables::A_T_RE_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return A_T_Re_binned(obs);
    case Observables::A_T_RE_CPV_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return A_T_Re_CPV_binned(obs);
    case Observables::A_IM_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return A_Im_binned(obs);
    case Observables::ALPHA_K_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return alpha_K_binned(obs);
    case Observables::H_T_1_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return H_T_1_binned(obs);
    case Observables::H_T_2_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return H_T_2_binned(obs);
    case Observables::H_T_3_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return H_T_3_binned(obs);
    case Observables::P_1_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return A_T_2_binned(obs);
    case Observables::P_2_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return P_2_binned(obs);
    case Observables::P_3_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);   
        return P_3_binned(obs);
    case Observables::P_4_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return H_T_1_binned(obs);
    case Observables::P_5_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return H_T_2_binned(obs);
    case Observables::P_6_B0__KSTAR0_TAU_TAU:   
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);   
        return P_6_binned(obs);
    case Observables::P_8_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return P_8_binned(obs);
    case Observables::P_PRIME_4_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(4, false, obs);
    case Observables::P_PRIME_5_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(5, false, obs);
    case Observables::P_PRIME_6_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(6, false, obs);
    case Observables::P_PRIME_8_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(8, false, obs);
    case Observables::S_1C_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(0, false, obs);
    case Observables::S_2S_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(1, false, obs);
    case Observables::S_2C_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(2, false, obs);
    case Observables::S_3_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(3, false, obs);
    case Observables::S_4_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(4, false, obs);
    case Observables::S_5_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(5, false, obs);
    case Observables::S_6C_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return S_i_binned(10, false, obs);
    case Observables::S_7_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(7, false, obs);
    case Observables::S_8_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(8, false, obs);
    case Observables::S_9_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(9, false, obs);
    case Observables::A_1C_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(0, true, obs);
    case Observables::A_2S_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(1, true, obs);
    case Observables::A_FL_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(2, true, obs);
    case Observables::A_3_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return S_i_binned(3, true, obs);
    case Observables::A_4_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(4, true, obs);
    case Observables::A_5_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(5, true, obs);
    case Observables::A_6S_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(6, true, obs);
    case Observables::A_6C_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(10, true, obs);
    case Observables::A_7_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(7, true, obs);
    case Observables::A_8_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(8, true, obs);
    case Observables::A_9_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return S_i_binned(9, true, obs);
    case Observables::P_1_CPV_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return P_i_CPV_binned(1, obs);
    case Observables::P_2_CPV_B0__KSTAR0_TAU_TAU:      
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);
        return P_i_CPV_binned(2, obs);
    case Observables::P_3_CPV_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return P_i_CPV_binned(3, obs);
    case Observables::P_PRIME_4_CPV_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(4, true, obs);
    case Observables::P_PRIME_5_CPV_B0__KSTAR0_TAU_TAU:    
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0);  
        return Pp_i_binned(5, true, obs);
    case Observables::P_PRIME_6_CPV_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(6, true, obs);
    case Observables::P_PRIME_8_CPV_B0__KSTAR0_TAU_TAU:     
        set_lepton_gen_and_charge(BKstarllConfig::Lepton::TAU, BKstarllConfig::B_Charge::B_0); 
        return Pp_i_binned(8, true, obs);
    case Observables::R_1_B__KSTAR_L_L:
        return Rm1_BKstar(obs, BKstarllConfig::B_Charge::B_PLUS);

    case Observables::R_1_B0__KSTAR0_L_L:
        return Rm1_BKstar(obs, BKstarllConfig::B_Charge::B_0);
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> BKstarllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}

bool BKstarllDecay::is_observable_binned(ObservableId obs) const {
    const auto enum_obs = ObservableMapper::enum_of(obs);
    if (!enum_obs.has_value()) {
        return true;
    }

    switch (enum_obs.value()) {
    case Observables::Q0_A_FB_B__KSTAR_E_E:
    case Observables::Q0_A_FB_B__KSTAR_MU_MU:
    case Observables::Q0_A_FB_B__KSTAR_TAU_TAU:
    case Observables::Q0_A_FB_B0__KSTAR0_E_E:
    case Observables::Q0_A_FB_B0__KSTAR0_MU_MU:
    case Observables::Q0_A_FB_B0__KSTAR0_TAU_TAU:
        return false;
    default:
        return true;
    }
}

void BKstarllDecay::set_n_threads(size_t n_threads) {
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