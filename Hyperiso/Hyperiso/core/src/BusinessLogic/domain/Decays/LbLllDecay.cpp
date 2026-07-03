#include "LbLllDecay.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>

void LbLllDecay::load_params() {
    fill_wilson_cache();

    cache.ff_calculator = LbLFFCalculator(p, cfg.ff_src);
    cache.m_b_mu_b = (*iobs_qcdp)(MassConfig(5, (*p)(ParamId{ParameterType::WILSON, "B_SCALE", 1}, DataType::VALUE), MassType::MSBAR, MassType::POLE));
    cache.m_Lb = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 5122}, DataType::VALUE);
    cache.m_L = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 3122}, DataType::VALUE);
    cache.life_L = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 5122}, DataType::VALUE) / HBAR;
    cache.alpha_L = (*p)(ParamId{ParameterType::DECAY, "Lb_L", 8}, DataType::VALUE);
    cache.N_0 = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE) * (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE) / (std::sqrt(6144. * std::pow(PI, 5) * std::pow(cache.m_Lb, 3)));
    cache.q2_max = std::pow(cache.m_Lb - cache.m_L, 2);

    load_cfg_dep_params();

    // printf("alpha_em = %.4e\n", (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE).real());
    // printf("N0 = %.4e + %.4e i\n", cache.N_0.real(), cache.N_0.imag());

    // printf("f_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::F_PERP, 1.0));
    // printf("h_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::H_PERP, 1.0));
    // printf("g_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::G_PERP, 1.0));
    // printf("f_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::F_PLUS, 1.0));
    // printf("h_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::H_PLUS, 1.0));
    // printf("g_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::G_PLUS, 1.0));
    // printf("h_tilde_perp = %.4e\n", cache.ff_calculator.get(LbL_FF::H_TILDE_PERP, 1.0));
    // printf("h_tilde_+ = %.4e\n", cache.ff_calculator.get(LbL_FF::H_TILDE_PLUS, 1.0));
}

void LbLllDecay::fill_wilson_cache() {
    cache.C.clear();

    auto b_wilsons  = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);

    cache.C[WCoef::C7]   = b_wilsons.at(WCoef::C7);
    cache.C[WCoef::C9]   = b_wilsons.at(WCoef::C9);
    cache.C[WCoef::C10]  = b_wilsons.at(WCoef::C10);

    cache.C[WCoef::CP7]  = bp_wilsons.at(WCoef::CP7);
    cache.C[WCoef::CP9]  = bp_wilsons.at(WCoef::CP9);
    cache.C[WCoef::CP10] = bp_wilsons.at(WCoef::CP10);
}

void LbLllDecay::set_cfg_flags(LbLllConfig::Lepton gen) {
    if (cfg.gen != gen) {
        cfg.gen = gen;
        load_cfg_dep_params();
    }
}

void LbLllDecay::load_cfg_dep_params() {
    cache.m_l = (*p)(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen}, DataType::VALUE);
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    compute_binned_K_i();
}

double LbLllDecay::beta_l(double q2) {
    return std::sqrt(1 - 4. * cache.m_l * cache.m_l / q2);
}

double LbLllDecay::lambda(double q2) {
    double mLb2 = cache.m_Lb * cache.m_Lb;
    double mL2 = cache.m_L * cache.m_L;
    return mLb2 * mLb2 + mL2 * mL2 + q2 * q2 - 2. * (mLb2 * mL2 + (mLb2 + mL2) * q2);
}

double LbLllDecay::s_p(double q2) {
    return std::pow(cache.m_Lb + cache.m_L, 2) - q2;
}

double LbLllDecay::s_m(double q2) {
    return std::pow(cache.m_Lb - cache.m_L, 2) - q2;
}

complex_t LbLllDecay::N(double q2, bool bar) {
    complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0; 
    return N0 * std::sqrt(q2 * beta_l(q2) * std::sqrt(lambda(q2)));
}

complex_t LbLllDecay::A_perp_1(double q2, double sign, bool bar) {
    complex_t C_V = (cache.C[WCoef::C9] + cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
    complex_t C_T = cache.C[WCoef::C7] + cache.C[WCoef::CP7];
    if (bar) {
        C_V = std::conj(C_V);
        C_T = std::conj(C_T);
    }
    complex_t HVplus = -cache.ff_calculator.get(LbL_FF::F_PERP, q2) * std::sqrt(2 * s_m(q2));
    complex_t HTplus = cache.ff_calculator.get(LbL_FF::H_PERP, q2) * (cache.m_Lb + cache.m_L) * std::sqrt(2 * s_m(q2));
    return RT2 * N(q2, bar) * (C_V * HVplus - 2 * cache.m_b_mu_b / q2 * C_T * HTplus);
}

complex_t LbLllDecay::A_par_1(double q2, double sign, bool bar) {
    complex_t C_A = (cache.C[WCoef::C9] - cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
    complex_t C_TA = cache.C[WCoef::C7] - cache.C[WCoef::CP7];
    if (bar) {
        C_A = std::conj(C_A);
        C_TA = std::conj(C_TA);
    }
    complex_t HAplus = -cache.ff_calculator.get(LbL_FF::G_PERP, q2) * std::sqrt(2 * s_p(q2));
    complex_t HT5plus = -cache.ff_calculator.get(LbL_FF::H_TILDE_PERP, q2) * (cache.m_Lb - cache.m_L) * std::sqrt(2 * s_p(q2));
    return -RT2 * N(q2, bar) * (C_A * HAplus + 2 * cache.m_b_mu_b / q2 * C_TA * HT5plus);
}

complex_t LbLllDecay::A_perp_0(double q2, double sign, bool bar) {
    complex_t C_V = (cache.C[WCoef::C9] + cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
    complex_t C_T = cache.C[WCoef::C7] + cache.C[WCoef::CP7];
    if (bar) {
        C_V = std::conj(C_V);
        C_T = std::conj(C_T);
    }
    complex_t HV0 = cache.ff_calculator.get(LbL_FF::F_PLUS, q2) * (cache.m_Lb + cache.m_L) * std::sqrt(s_m(q2) / q2);
    complex_t HT0 = -cache.ff_calculator.get(LbL_FF::H_PLUS, q2) * std::sqrt(q2 * s_m(q2));
    return RT2 * N(q2, bar) * (C_V * HV0 - 2 * cache.m_b_mu_b / q2 * C_T * HT0);
}

complex_t LbLllDecay::A_par_0(double q2, double sign, bool bar) {
    complex_t C_A = (cache.C[WCoef::C9] - cache.C[WCoef::CP9]) + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
    complex_t C_TA = cache.C[WCoef::C7] - cache.C[WCoef::CP7];
    if (bar) {
        C_A = std::conj(C_A);
        C_TA = std::conj(C_TA);
    }
    complex_t HA0 = cache.ff_calculator.get(LbL_FF::G_PLUS, q2) * (cache.m_Lb - cache.m_L) * std::sqrt(s_p(q2) / q2);
    complex_t HT50 = cache.ff_calculator.get(LbL_FF::H_TILDE_PLUS, q2) * std::sqrt(q2 * s_p(q2));
    return -RT2 * N(q2, bar) * (C_A * HA0 + 2 * cache.m_b_mu_b / q2 * C_TA * HT50);
}

double LbLllDecay::K1ss(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return std::real(0.25*(ALpar1*conj(ALpar1)+ALperp1*conj(ALperp1)+ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1))
         + 0.25*(1.+b_l*b_l)*(ALpar0*conj(ALpar0)+ALperp0*conj(ALperp0)+ARpar0*conj(ARpar0)+ARperp0*conj(ARperp0))
         + 0.5*(1.-b_l*b_l)*(ARpar1*conj(ALpar1)+ARperp1*conj(ALperp1)+ARpar0*conj(ALpar0)+ARperp0*conj(ALperp0)));
}

double LbLllDecay::K1cc(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return std::real(0.25*(1.+b_l*b_l)*(ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1)+ALpar1*conj(ALpar1)+ALperp1*conj(ALperp1))
					+0.25*(1.-b_l*b_l)*(ARpar0*conj(ARpar0)+ARperp0*conj(ARperp0)+ALpar0*conj(ALpar0)+ALperp0*conj(ALperp0))
					+0.5*(1.-b_l*b_l)*(ARpar1*conj(ALpar1)+ARperp1*conj(ALperp1)+ARpar0*conj(ALpar0)+ARperp0*conj(ALperp0)));
}

double LbLllDecay::K1c(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    double b_l = beta_l(q2);
    return -b_l*std::real(ARperp1*conj(ARpar1)-ALperp1*conj(ALpar1));;
}

double LbLllDecay::K2ss(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return 0.5 * cache.alpha_L * std::real(
        ARperp1*conj(ARpar1)+ALperp1*conj(ALpar1)
      + (1. + b_l*b_l)*(ARperp0*conj(ARpar0)+ALperp0*conj(ALpar0))
      + (1. - b_l*b_l)*(ARperp1*conj(ALpar1)+ARpar1*conj(ALperp1)+ARperp0*conj(ALpar0)+ARpar0*conj(ALperp0)));
}

double LbLllDecay::K2cc(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    complex_t ALpar0 = A_par_0(q2, -1, bar);
    complex_t ALperp0 = A_perp_0(q2, -1, bar);
    complex_t ARpar0 = A_par_0(q2, 1, bar);
    complex_t ARperp0 = A_perp_0(q2, 1, bar);

    double b_l = beta_l(q2);
    
    return 0.5 * cache.alpha_L * std::real(
        (1.+b_l*b_l)*(ARperp1*conj(ARpar1)+ALperp1*conj(ALpar1))
      + (1.-b_l*b_l)*(ARpar0*conj(ARperp0)+ALpar0*conj(ALperp0))
      + (1.-b_l*b_l)*(ARperp1*conj(ALpar1)+ARpar1*conj(ALperp1)+ARperp0*conj(ALpar0)+ARpar0*conj(ALperp0)));
}

double LbLllDecay::K2c(double q2, bool bar) {
    complex_t ALpar1 = A_par_1(q2, -1, bar);
    complex_t ALperp1 = A_perp_1(q2, -1, bar);
    complex_t ARpar1 = A_par_1(q2, 1, bar);
    complex_t ARperp1 = A_perp_1(q2, 1, bar);
    double b_l = beta_l(q2);
    return -0.5*cache.alpha_L*b_l*std::real(ARpar1*conj(ARpar1)+ARperp1*conj(ARperp1)-ALpar1*conj(ALpar1)-ALperp1*conj(ALperp1));
}

void LbLllDecay::compute_binned_K_i() {
    // Same strategy as BKstarllDecay::compute_binned_J_i:
    // - one quadrature pass per bin and CP-conjugation
    // - evaluate the 8 transversity amplitudes once per q2 point
    // - derive all K_i from those amplitudes
    // - split independent bins across cfg.n_threads workers
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
        complex_t ALpar1;
        complex_t ALperp1;
        complex_t ARpar1;
        complex_t ARperp1;
        complex_t ALpar0;
        complex_t ALperp0;
        complex_t ARpar0;
        complex_t ARperp0;
        double beta;
        double beta2;
    };

    const auto& bins = this->bins.value();
    const size_t nbins = bins.size();

    auto clear_and_resize = [&] (std::array<std::vector<double>, 6>& dest) {
        for (auto& v : dest) {
            v.assign(nbins, std::numeric_limits<double>::quiet_NaN());
        }
    };

    clear_and_resize(cache.K_i_binned);
    clear_and_resize(cache.K_i_bar_binned);
    cache.bin_widths.assign(nbins, std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < nbins; ++i) {
        cache.bin_widths[i] = bins[i].second - bins[i].first;
    }

    if (nbins == 0) {
        return;
    }

    auto eval_amplitudes = [&] (double q2, bool bar) -> AmpSet {
        const double beta = beta_l(q2);
        return {
            A_par_1(q2,  -1, bar),
            A_perp_1(q2, -1, bar),
            A_par_1(q2,   1, bar),
            A_perp_1(q2,  1, bar),
            A_par_0(q2,  -1, bar),
            A_perp_0(q2, -1, bar),
            A_par_0(q2,   1, bar),
            A_perp_0(q2,  1, bar),
            beta,
            beta * beta
        };
    };

    auto eval_integrands = [&] (double q2, bool bar) -> std::array<double, 6> {
        const AmpSet a = eval_amplitudes(q2, bar);

        const double norm_ALpar1  = std::norm(a.ALpar1);
        const double norm_ALperp1 = std::norm(a.ALperp1);
        const double norm_ARpar1  = std::norm(a.ARpar1);
        const double norm_ARperp1 = std::norm(a.ARperp1);
        const double norm_ALpar0  = std::norm(a.ALpar0);
        const double norm_ALperp0 = std::norm(a.ALperp0);
        const double norm_ARpar0  = std::norm(a.ARpar0);
        const double norm_ARperp0 = std::norm(a.ARperp0);

        const double one_plus_beta2 = 1.0 + a.beta2;
        const double one_minus_beta2 = 1.0 - a.beta2;
        const double alpha_L = cache.alpha_L;

        const complex_t k1_cross =
            a.ARpar1  * std::conj(a.ALpar1)
          + a.ARperp1 * std::conj(a.ALperp1)
          + a.ARpar0  * std::conj(a.ALpar0)
          + a.ARperp0 * std::conj(a.ALperp0);

        const complex_t k2_cross =
            a.ARperp1 * std::conj(a.ALpar1)
          + a.ARpar1  * std::conj(a.ALperp1)
          + a.ARperp0 * std::conj(a.ALpar0)
          + a.ARpar0  * std::conj(a.ALperp0);

        const double k1ss = std::real(
            0.25 * (norm_ALpar1 + norm_ALperp1 + norm_ARpar1 + norm_ARperp1)
          + 0.25 * one_plus_beta2 * (norm_ALpar0 + norm_ALperp0 + norm_ARpar0 + norm_ARperp0)
          + 0.5  * one_minus_beta2 * k1_cross
        );

        const double k1cc = std::real(
            0.25 * one_plus_beta2 * (norm_ARpar1 + norm_ARperp1 + norm_ALpar1 + norm_ALperp1)
          + 0.25 * one_minus_beta2 * (norm_ARpar0 + norm_ARperp0 + norm_ALpar0 + norm_ALperp0)
          + 0.5  * one_minus_beta2 * k1_cross
        );

        const double k1c = -a.beta * std::real(
            a.ARperp1 * std::conj(a.ARpar1)
          - a.ALperp1 * std::conj(a.ALpar1)
        );

        const double k2ss = 0.5 * alpha_L * std::real(
            a.ARperp1 * std::conj(a.ARpar1)
          + a.ALperp1 * std::conj(a.ALpar1)
          + one_plus_beta2 * (a.ARperp0 * std::conj(a.ARpar0) + a.ALperp0 * std::conj(a.ALpar0))
          + one_minus_beta2 * k2_cross
        );

        const double k2cc = 0.5 * alpha_L * std::real(
            one_plus_beta2 * (a.ARperp1 * std::conj(a.ARpar1) + a.ALperp1 * std::conj(a.ALpar1))
          + one_minus_beta2 * (a.ARpar0 * std::conj(a.ARperp0) + a.ALpar0 * std::conj(a.ALperp0))
          + one_minus_beta2 * k2_cross
        );

        const double k2c = -0.5 * alpha_L * a.beta * (
            norm_ARpar1 + norm_ARperp1 - norm_ALpar1 - norm_ALperp1
        );

        return {{k1ss, k1cc, k1c, k2ss, k2cc, k2c}};
    };

    auto integrate_bin = [&] (double q2_l, double q2_u, bool bar) -> std::array<double, 6> {
        std::array<double, 6> acc {};

        const double width = q2_u - q2_l;
        if (!std::isfinite(width) || width <= 0.0) {
            acc.fill(std::numeric_limits<double>::quiet_NaN());
            return acc;
        }

        const double center = 0.5 * (q2_l + q2_u);
        const double half_width = 0.5 * width;

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

    auto compute_and_store_bin = [&] (size_t bin_idx) {
        const auto& [q2_l, q2_u] = bins[bin_idx];

        const auto integ = integrate_bin(q2_l, q2_u, false);
        const auto integ_bar = integrate_bin(q2_l, q2_u, true);

        for (size_t k = 0; k < cache.K_i_binned.size(); ++k) {
            cache.K_i_binned[k][bin_idx] = integ[k];
            cache.K_i_bar_binned[k][bin_idx] = integ_bar[k];
        }
    };

    auto requested_threads = cfg.n_threads;
    if (requested_threads == 0u) {
        requested_threads = std::thread::hardware_concurrency();
    }
    if (requested_threads == 0u) {
        requested_threads = 1u;
    }

    const size_t nworkers = std::min<size_t>(requested_threads, nbins);

    if (nworkers <= 1u) {
        for (size_t i = 0; i < nbins; ++i) {
            compute_and_store_bin(i);
        }
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(nworkers);

    std::exception_ptr first_exception = nullptr;
    std::mutex exception_mutex;

    auto worker = [&] (size_t begin, size_t end) {
        try {
            for (size_t i = begin; i < end; ++i) {
                compute_and_store_bin(i);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!first_exception) {
                first_exception = std::current_exception();
            }
        }
    };

    const size_t chunk = (nbins + nworkers - 1) / nworkers;
    for (size_t w = 0; w < nworkers; ++w) {
        const size_t begin = w * chunk;
        const size_t end = std::min(nbins, begin + chunk);
        if (begin >= end) {
            break;
        }
        workers.emplace_back(worker, begin, end);
    }

    for (auto& th : workers) {
        th.join();
    }

    if (first_exception) {
        std::rethrow_exception(first_exception);
    }
}

std::vector<ObservableValue> LbLllDecay::dBR_dq2_binned(Observables oid, bool br) {
    std::vector<ObservableValue> out;
    double br_factor = br ? cache.life_L : 1.0;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1ss = cache.K_i_binned[0][i] + cache.K_i_bar_binned[0][i];
        double K1cc = cache.K_i_binned[1][i] + cache.K_i_bar_binned[1][i];
        double integrated_rate = (2 * K1ss + K1cc) / 2;
        const double requested_width = this->bins.value()[i].second - this->bins.value()[i].first;
        const double width =
            (i < cache.bin_widths.size() && std::isfinite(cache.bin_widths[i]) && cache.bin_widths[i] > 0.0)
            ? cache.bin_widths[i]
            : requested_width;

        const double res =
            (std::isfinite(width) && width > 0.0)
            ? integrated_rate * br_factor / width
            : std::numeric_limits<double>::quiet_NaN();

        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }

    return out;
}

double LbLllDecay::dG_dq2_avg_bin(size_t bin) {
    double K1ss = cache.K_i_binned[0][bin] + cache.K_i_bar_binned[0][bin];
    double K1cc = cache.K_i_binned[1][bin] + cache.K_i_bar_binned[1][bin];
    return (2 * K1ss + K1cc) / 2;
}

std::vector<ObservableValue> LbLllDecay::A_FB_l(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1c = cache.K_i_binned[2][i] + cache.K_i_bar_binned[2][i];
        double res = 0.75 * K1c / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::A_FB_h(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K2ss = cache.K_i_binned[3][i] + cache.K_i_bar_binned[3][i];
        double K2cc = cache.K_i_binned[4][i] + cache.K_i_bar_binned[4][i];
        double res = 0.25 * (2 * K2ss + K2cc) / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::A_FB_lh(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K2c = cache.K_i_binned[5][i] + cache.K_i_bar_binned[5][i];
        double res = 0.375 * K2c / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::F_L(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1ss = cache.K_i_binned[0][i] + cache.K_i_bar_binned[0][i];
        double K1cc = cache.K_i_binned[1][i] + cache.K_i_bar_binned[1][i];
        double res = 0.5 * (2 * K1ss - K1cc) / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::F_T(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double K1cc = cache.K_i_binned[1][i] + cache.K_i_bar_binned[1][i];
        double res = K1cc / dG_dq2_avg_bin(i);
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> LbLllDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::DBR_DQ2_LAMBDA_B__LAMBDA_E_E:   
        set_cfg_flags(LbLllConfig::Lepton::E);
        return dBR_dq2_binned(obs, true);
    case Observables::DGAMMA_DQ2_LAMBDA_B__LAMBDA_E_E:   
        set_cfg_flags(LbLllConfig::Lepton::E);
        return dBR_dq2_binned(obs, false);
    case Observables::A_FB_L_LAMBDA_B__LAMBDA_E_E:   
        set_cfg_flags(LbLllConfig::Lepton::E);
        return A_FB_l(obs);
    case Observables::A_FB_H_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return A_FB_h(obs);
    case Observables::A_FB_LH_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return A_FB_lh(obs);
    case Observables::F_L_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return F_L(obs);
    case Observables::F_T_LAMBDA_B__LAMBDA_E_E:
        set_cfg_flags(LbLllConfig::Lepton::E);
        return F_T(obs);
    case Observables::DBR_DQ2_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);
        return dBR_dq2_binned(obs, true);
    case Observables::DGAMMA_DQ2_LAMBDA_B__LAMBDA_MU_MU:   
        set_cfg_flags(LbLllConfig::Lepton::MU);
        return dBR_dq2_binned(obs, false);
    case Observables::A_FB_L_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);
        return A_FB_l(obs);
    case Observables::A_FB_H_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return A_FB_h(obs);
    case Observables::A_FB_LH_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return A_FB_lh(obs);
    case Observables::F_L_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return F_L(obs);
    case Observables::F_T_LAMBDA_B__LAMBDA_MU_MU:
        set_cfg_flags(LbLllConfig::Lepton::MU);   
        return F_T(obs);
    case Observables::DBR_DQ2_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return dBR_dq2_binned(obs, true);
    case Observables::DGAMMA_DQ2_LAMBDA_B__LAMBDA_TAU_TAU:   
        set_cfg_flags(LbLllConfig::Lepton::TAU);
        return dBR_dq2_binned(obs, false);
    case Observables::A_FB_L_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return A_FB_l(obs);
    case Observables::A_FB_H_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return A_FB_h(obs);
    case Observables::A_FB_LH_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return A_FB_lh(obs);
    case Observables::F_L_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return F_L(obs);
    case Observables::F_T_LAMBDA_B__LAMBDA_TAU_TAU:
        set_cfg_flags(LbLllConfig::Lepton::TAU);   
        return F_T(obs);
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> LbLllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}


void LbLllDecay::set_n_threads(size_t n_threads) {
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
