#include "BKllDecay.h"

#include <algorithm>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

using Charge = BKllConfig::B_Charge;

const std::unordered_set<ObservableId> BKllDecay::dBR_dq2_ids = {
    ObservableMapper::to_id(Observables::DBR_DQ2_B0__K0_E_E), 
    ObservableMapper::to_id(Observables::DBR_DQ2_B0__K0_MU_MU), 
    ObservableMapper::to_id(Observables::DBR_DQ2_B0__K0_TAU_TAU), 
    ObservableMapper::to_id(Observables::DBR_DQ2_B__K_E_E), 
    ObservableMapper::to_id(Observables::DBR_DQ2_B__K_MU_MU), 
    ObservableMapper::to_id(Observables::DBR_DQ2_B__K_TAU_TAU)
};

const std::unordered_set<ObservableId> BKllDecay::dG_dq2_ids = {
    ObservableMapper::to_id(Observables::DGAMMA_DQ2_B0__K0_E_E), 
    ObservableMapper::to_id(Observables::DGAMMA_DQ2_B0__K0_MU_MU), 
    ObservableMapper::to_id(Observables::DGAMMA_DQ2_B0__K0_TAU_TAU), 
    ObservableMapper::to_id(Observables::DGAMMA_DQ2_B__K_E_E), 
    ObservableMapper::to_id(Observables::DGAMMA_DQ2_B__K_MU_MU), 
    ObservableMapper::to_id(Observables::DGAMMA_DQ2_B__K_TAU_TAU)
};

const std::unordered_set<ObservableId> BKllDecay::A_FB_ids = {
    ObservableMapper::to_id(Observables::A_FB_B0__K0_E_E), 
    ObservableMapper::to_id(Observables::A_FB_B0__K0_MU_MU), 
    ObservableMapper::to_id(Observables::A_FB_B0__K0_TAU_TAU), 
    ObservableMapper::to_id(Observables::A_FB_B__K_E_E), 
    ObservableMapper::to_id(Observables::A_FB_B__K_MU_MU), 
    ObservableMapper::to_id(Observables::A_FB_B__K_TAU_TAU)
};

const std::unordered_set<ObservableId> BKllDecay::F_H_ids = {
    ObservableMapper::to_id(Observables::F_H_B0__K0_E_E), 
    ObservableMapper::to_id(Observables::F_H_B0__K0_MU_MU), 
    ObservableMapper::to_id(Observables::F_H_B0__K0_TAU_TAU), 
    ObservableMapper::to_id(Observables::F_H_B__K_E_E), 
    ObservableMapper::to_id(Observables::F_H_B__K_MU_MU), 
    ObservableMapper::to_id(Observables::F_H_B__K_TAU_TAU)
};


const std::map<Observables, std::pair<BKllConfig::Lepton, BKllConfig::B_Charge>> BKllDecay::cfg_map {
    {Observables::DGAMMA_DQ2_B0__K0_E_E, {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_0}},
    {Observables::DBR_DQ2_B0__K0_E_E,    {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_0}},
    {Observables::A_FB_B0__K0_E_E,       {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_0}},
    {Observables::F_H_B0__K0_E_E,        {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_0}},

    {Observables::DGAMMA_DQ2_B__K_E_E, {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_PLUS}},
    {Observables::DBR_DQ2_B__K_E_E,    {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_PLUS}},
    {Observables::A_FB_B__K_E_E,       {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_PLUS}},
    {Observables::F_H_B__K_E_E,        {BKllConfig::Lepton::E, BKllConfig::B_Charge::B_PLUS}},

    {Observables::DGAMMA_DQ2_B0__K0_MU_MU, {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_0}},
    {Observables::DBR_DQ2_B0__K0_MU_MU,    {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_0}},
    {Observables::A_FB_B0__K0_MU_MU,       {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_0}},
    {Observables::F_H_B0__K0_MU_MU,        {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_0}},

    {Observables::DGAMMA_DQ2_B__K_MU_MU, {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_PLUS}},
    {Observables::DBR_DQ2_B__K_MU_MU,    {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_PLUS}},
    {Observables::A_FB_B__K_MU_MU,       {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_PLUS}},
    {Observables::F_H_B__K_MU_MU,        {BKllConfig::Lepton::MU, BKllConfig::B_Charge::B_PLUS}},

    {Observables::DGAMMA_DQ2_B0__K0_TAU_TAU, {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_0}},
    {Observables::DBR_DQ2_B0__K0_TAU_TAU,    {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_0}},
    {Observables::A_FB_B0__K0_TAU_TAU,       {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_0}},
    {Observables::F_H_B0__K0_TAU_TAU,        {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_0}},

    {Observables::DGAMMA_DQ2_B__K_TAU_TAU, {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_PLUS}},
    {Observables::DBR_DQ2_B__K_TAU_TAU,    {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_PLUS}},
    {Observables::A_FB_B__K_TAU_TAU,       {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_PLUS}},
    {Observables::F_H_B__K_TAU_TAU,        {BKllConfig::Lepton::TAU, BKllConfig::B_Charge::B_PLUS}},
};

void BKllDecay::load_params() {
    fill_wilson_cache();

    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.G_F = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    cache.m_s = (*p)(ParamId{ParameterType::SM, "MASS", 3}, DataType::VALUE);
    cache.mu_b = (*p)(ParamId{ParameterType::WILSON, "B_SCALE", 1}, DataType::VALUE);
    cache.alpha_s_mu_b = (*iobs_qcdp)(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    // cache.m_c_mu_b = (*iobs_qcdp)(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_c_mu_b = (*p)(ParamId{ParameterType::SM, "MASS", 4}, DataType::VALUE); // ASK: To match Superiso
    cache.m_b_mu_b = (*iobs_qcdp)(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_b_m_b = (*p)(ParamId{ParameterType::SM, "QCD", {5, 1}}, DataType::VALUE); // ASK: To match SI at high q² : why not m_b(mu_b) ?
    double mu_f = sqrt(cache.mu_b * (*p)(ParamId{ParameterType::DECAY, "B_K", 14}, DataType::VALUE));
    cache.m_b_PS = (*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE) - 4 * (*iobs_qcdp)(AlphasConfig((*p)(ParamId{ParameterType::SM, "QCD", {5, 2}}, DataType::VALUE), MassType::POLE, MassType::POLE)) * mu_f / (3 * PI);
    cache.L_b = std::log(cache.mu_b / cache.m_b_PS);
    cache.Delta_M = -6. * cache.L_b - 4. * (1 - mu_f / cache.m_b_PS);
    cache.lambda_hat_u = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE) 
                            / (std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE));
    cache.q2_low = (*p)(ParamId{ParameterType::DECAY, "B_K", {15, 1}}, DataType::VALUE);
    cache.q2_high = (*p)(ParamId{ParameterType::DECAY, "B_K", {15, 2}}, DataType::VALUE);

    for (size_t i = 0; i < 4; i++) {
        cache.A_had_err_low_0[i] = (*p)(ParamId{ParameterType::DECAY, "B_K", {18, 1, i + 1}}, DataType::VALUE);
        cache.A_had_err_low_1[i] = (*p)(ParamId{ParameterType::DECAY, "B_K", {18, 2, i + 1}}, DataType::VALUE);
        cache.A_had_err_high[i] = (*p)(ParamId{ParameterType::DECAY, "B_K", {18, 3, i + 1}}, DataType::VALUE);
    }

    load_cfg_dependent_params();

    // printf("alpha_em = %.4e\n", cache.alpha_em);
    // printf("m_l = %.4e\n", cache.m_l);
    // printf("m_s = %.4e\n", cache.m_s);
    // printf("mu_b = %.4e\n", cache.mu_b);
    // printf("alpha_s(mu_b) = %.4e\n", cache.alpha_s_mu_b);
    // printf("m_c(mu_b) = %.4e\n", cache.m_c_mu_b);
    // printf("m_b(mu_b) = %.4e\n", cache.m_b_mu_b);
    // printf("m_b_PS = %.4e\n", cache.m_b_PS);
    // printf("L_b = %.4e\n", cache.L_b);
    // printf("m_B = %.4e\n", cache.m_B);
    // printf("m_K = %.4e\n", cache.m_K);
    // printf("Delta_M = %.4e\n", cache.Delta_M);
    // printf("lambda_hat_u = %.4e + %.4e i\n", cache.lambda_hat_u.real(), cache.lambda_hat_u.imag());
    // printf("N_0 = %.4e\n", cache.N_0);

    // printf("f_0(s = 1.0 GeV²) = %.4e\n", cache.ff_calculator.get(BP_FF::F_0, 1.0));
    // printf("f_+(s = 1.0 GeV²) = %.4e\n", cache.ff_calculator.get(BP_FF::F_PLUS, 1.0));
    // printf("f_T(s = 1.0 GeV²) = %.4e\n", cache.ff_calculator.get(BP_FF::F_T, 1.0));

    // printf("T_P = %.4e + %.4e i\n", real(cache.qcdf_calculator.T_P(1.0, false)), imag(cache.qcdf_calculator.T_P(1.0, false)));

    // printf("F_A(s = 1.0 GeV²) = %.4e + %.4e i\n", real(F_A(1.0)), imag(F_A(1.0)));
    // printf("F_V(s = 1.0 GeV²) = %.4e + %.4e i\n", real(F_V(1.0)), imag(F_V(1.0)));
    // printf("F_S(s = 1.0 GeV²) = %.4e + %.4e i\n", real(F_S(1.0)), imag(F_S(1.0)));
    // printf("F_P(s = 1.0 GeV²) = %.4e + %.4e i\n", real(F_P(1.0)), imag(F_P(1.0)));
}

void BKllDecay::fill_wilson_cache() {
    cache.C.clear();

    auto b_wilsons  = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);

    WCoef bp_cached[6] {
        WCoef::CP7, WCoef::CP8, WCoef::CP9,
        WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2
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

void BKllDecay::load_cfg_dependent_params() {
    const int B_id = cfg.charge == Charge::B_0 ? 511 : 521;
    const int P_id = cfg.charge == Charge::B_0 ? 311 : 321;

    cache.ff_calculator = BPFFCalculator(
        B_id,
        P_id,
        p,
        cfg.ff_src
    );

    cache.qcdf_calculator = BPQCDfCalculator(
        B_id,
        P_id,
        cache.mu_b,
        cache.C,
        std::make_shared<BPFFCalculator>(cache.ff_calculator),
        cfg.ff_type,
        p,
        iobs_qcdp
    );

    cache.m_l = (*p)(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen}, DataType::VALUE);
    cache.m_B = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", B_id}, DataType::VALUE);
    cache.m_K = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", P_id}, DataType::VALUE);
    cache.life_B = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", B_id}, DataType::VALUE) / HBAR;
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    cache.q2_max = std::pow(cache.m_B - cache.m_K, 2);
    cache.q2_lookup_min = std::max(cache.q2_min, 1e-4);
    cache.N_0 = std::pow(std::abs(std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE)) * cache.G_F * cache.alpha_em, 2) / (512. * std::pow(PI, 5) * std::pow(cache.m_B, 3));

    auto requested_threads = cfg.n_threads;
    if (requested_threads == 0u) {
        requested_threads = std::thread::hardware_concurrency();
    }
    if (requested_threads == 0u) {
        requested_threads = 1u;
    }

    const size_t npts = BKllCache::LOOKUP_SIZE;
    const size_t nworkers = std::min<size_t>(requested_threads, npts);

    if (nworkers <= 1u) {
        auto lam_T_P = [this] (double q2, bool bar) {
            return cache.qcdf_calculator.T_P(q2, bar);
        };

        fill_cache(lam_T_P, cache.q2_lookup_min, cache.q2_high, cache.T_P_lookup, false);
    } else {
        const double x_min = cache.q2_lookup_min;
        const double x_max = cache.q2_high;
        const double step = (x_max - x_min) / static_cast<double>(npts - 1);

        std::vector<std::shared_ptr<BPQCDfCalculator>> qcdf_locals;
        qcdf_locals.reserve(nworkers);
        for (size_t w = 0; w < nworkers; ++w) {
            auto ff_local = std::make_shared<BPFFCalculator>(cache.ff_calculator);
            qcdf_locals.emplace_back(std::make_shared<BPQCDfCalculator>(
                B_id,
                P_id,
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
                BPQCDfCalculator& qcdf = *qcdf_locals[worker_id];
                for (size_t i = begin; i < end; ++i) {
                    const double q2 = x_min + step * static_cast<double>(i);
                    cache.T_P_lookup[i] = qcdf.T_P(q2, false);
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

    compute_binned_abc();
}

void BKllDecay::set_lepton_gen_and_charge(BKllConfig::Lepton gen, BKllConfig::B_Charge charge) {
    bool changed = cfg.gen != gen || cfg.charge != charge;
    if (changed) {
        cfg.gen = gen;  
        cfg.charge = charge;
        load_cfg_dependent_params();
    }
}

complex_t BKllDecay::T_P_cached(double q2) {
    const double x = std::max(q2, cache.q2_lookup_min);
    return lerp(x, cache.T_P_lookup, cache.q2_lookup_min, cache.q2_high);
}

double BKllDecay::beta_l(double q2) {
    const double x = 1.0 - std::pow(2.0 * cache.m_l, 2) / q2;
    return std::sqrt(std::max(0.0, x));
}

double BKllDecay::lambda(double q2) {
    const double mB2 = cache.m_B * cache.m_B;
    const double mK2 = cache.m_K * cache.m_K;

    const double lam =
        mB2 * mB2
        + mK2 * mK2
        + q2 * q2
        - 2.0 * (mB2 * mK2 + (mB2 + mK2) * q2);

    return std::max(0.0, lam);
}

double BKllDecay::N(double q2) {
    return cache.N_0 * std::sqrt(lambda(q2)) * beta_l(q2);
}

complex_t BKllDecay::F_V_low(double q2) {
    complex_t F, F_T;
    double m_b_local;
    complex_t had_err_factor {1.0};

    if (cfg.ff_type == B_FF_Type::SOFT) {
        had_err_factor = 1.0 + cache.A_had_err_low_0[0] + cache.A_had_err_low_1[0] * q2 / 6.0;
        F_T = T_P_cached(q2);
        F = (cache.C[WCoef::C9] + cache.C[WCoef::CP9]) * cache.ff_calculator.get(BP_FF::XI_P, q2);
        m_b_local = cache.m_b_PS;
    } else {
        complex_t T_p_err_factor = 1.0 + cache.A_had_err_low_0[0] + cache.A_had_err_low_1[0] * q2 / 6.0;
        F_T = (cache.C[WCoef::C7] + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BP_FF::F_T, q2) + T_P_cached(q2) * T_p_err_factor;
        F = (cache.C[WCoef::C9] + cache.C[WCoef::CP9] + cache.qcdf_calculator.Y(q2)) * cache.ff_calculator.get(BP_FF::F_PLUS, q2);
        m_b_local = cache.m_b_PS + cache.alpha_s_mu_b * cache.Delta_M / (3 * PI);
    }

    // printf("C7 = %.4e\n", real(cache.C[WCoef::C7]));
    // printf("C9 = %.4e\n", real(cache.C[WCoef::C9]));
    // printf("C'7 = %.4e\n", real(cache.C[WCoef::CP7]));
    // printf("C'9 = %.4e\n", real(cache.C[WCoef::CP9]));
    // printf("m_b_local = %.4e\n", m_b_local);
    // printf("F = %.4e + %.4e i\n", real(F), imag(F));
    // printf("F_T = %.4e + %.4e i\n", real(F_T), imag(F_T));

    return (F + 2. * m_b_local / (cache.m_B + cache.m_K) * F_T) * had_err_factor;
}

complex_t BKllDecay::F_A_low(double q2) {
    double ff;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        ff = cache.ff_calculator.get(BP_FF::XI_P, q2);
    } else {
        ff = cache.ff_calculator.get(BP_FF::F_PLUS, q2);
    }

    return (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * ff;
}

complex_t BKllDecay::F_P_low(double q2) {
    double f_p, f0_fp;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        f_p = cache.ff_calculator.get(BP_FF::XI_P, q2);
        f0_fp = (std::pow(cache.m_B, 2) + std::pow(cache.m_K, 2) - q2) / std::pow(cache.m_B, 2) * cache.qcdf_calculator.Delta_P_0(q2);
    } else {
        f_p = cache.ff_calculator.get(BP_FF::F_PLUS, q2);
        f0_fp = cache.ff_calculator.get(BP_FF::F_0, q2) / f_p;
    }

    // printf("Delta_P_0(s = %.5f GeV²) = %.4e\n", q2, cache.qcdf_calculator.Delta_P_0(q2));
	// printf("f0_fp(s = %.5f GeV²) = %.4e\n", q2, f0_fp);

    return f_p * ((std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_mu_b - cache.m_s)) * f0_fp * (cache.C[WCoef::CQ2] + cache.C[WCoef::CPQ2])
                    - cache.m_l * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * (1. - (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / q2 * (f0_fp - 1.)));
}

complex_t BKllDecay::F_S_low(double q2) {
    double ff;

    if (cfg.ff_type == B_FF_Type::SOFT) {
        double f0_fp = 2. * (std::pow(cache.m_B, 2) + std::pow(cache.m_K, 2) - q2) / (2 * std::pow(cache.m_B, 2)) * cache.qcdf_calculator.Delta_P_0(q2);
        ff = cache.ff_calculator.get(BP_FF::XI_P, q2) * f0_fp;
    } else {
        ff = cache.ff_calculator.get(BP_FF::F_0, q2);
    }

    return (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_mu_b - cache.m_s)) * (cache.C[WCoef::CQ1] + cache.C[WCoef::CPQ1]) * ff;
}

complex_t BKllDecay::C7_eff(double q2) {
    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t A = BV::A_Seidel(s_hat, cache.L_b);
    return cache.C[WCoef::C7] + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*BV::f_87(s_hat, cache.L_b));
}

complex_t BKllDecay::C9_eff(double q2) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+cache.lambda_hat_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    double s_hat = q2 / std::pow(cache.m_b_PS, 2);
    complex_t B = BV::B_Seidel(s_hat, cache.L_b);
    complex_t C = BV::C_Seidel(q2, cache.mu_b);

    return cache.C[WCoef::C9]
         + BV::h(q2, 0., cache.mu_b) * C_h0
         + BV::h(q2, cache.m_b_PS, cache.mu_b) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * BV::f_89(s_hat))
         + std::pow(cache.m_c_mu_b, 2) / q2 * C_mc;
}

complex_t BKllDecay::F_V_high(double q2) {
    return (C9_eff(q2) + cache.C[WCoef::CP9]) * cache.ff_calculator.get(BP_FF::F_PLUS, q2) 
            + 2. * cache.m_b_m_b / (cache.m_B + cache.m_K) * (C7_eff(q2) + cache.C[WCoef::CP7]) * cache.ff_calculator.get(BP_FF::F_T, q2);
}

complex_t BKllDecay::F_A_high(double q2) {
    return (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * cache.ff_calculator.get(BP_FF::F_PLUS, q2);
}

complex_t BKllDecay::F_P_high(double q2) {
    double f_p = cache.ff_calculator.get(BP_FF::F_PLUS, q2);
    double f0_fp = cache.ff_calculator.get(BP_FF::F_0, q2) / f_p;
    return f_p * ((std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_m_b - cache.m_s)) * f0_fp * (cache.C[WCoef::CQ2] + cache.C[WCoef::CPQ2])
                    - cache.m_l * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]) * (1. - (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / q2 * (f0_fp - 1.)));
}

complex_t BKllDecay::F_S_high(double q2) {
    double ff = cache.ff_calculator.get(BP_FF::F_0, q2);
    return (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2)) / (2 * (cache.m_b_m_b - cache.m_s)) * (cache.C[WCoef::CQ1] + cache.C[WCoef::CPQ1]) * ff;
}

complex_t BKllDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BKllDecay::F_V(double q2) {
    return interpolate(q2, F_V_low(q2), F_V_high(q2));
}

complex_t BKllDecay::F_A(double q2) {
    return interpolate(q2, F_A_low(q2), F_A_high(q2));
}

complex_t BKllDecay::F_P(double q2) {
    return interpolate(q2, F_P_low(q2), F_P_high(q2));
}

complex_t BKllDecay::F_S(double q2) {
    return interpolate(q2, F_S_low(q2), F_S_high(q2));
}

double BKllDecay::a(double q2) {
    return N(q2) * (
        q2 * (std::pow(beta_l(q2) * std::abs(F_S(q2)), 2) + std::pow(std::abs(F_P(q2)), 2))
      + 0.25 * lambda(q2) * (std::pow(std::abs(F_A(q2)), 2) + std::pow(std::abs(F_V(q2)), 2))
      + 2 * cache.m_l * (std::pow(cache.m_B, 2) - std::pow(cache.m_K, 2) + q2) * std::real(F_P(q2) * std::conj(F_A(q2)))
      + std::pow(2 * cache.m_l * cache.m_B, 2) * std::pow(std::abs(F_A(q2)), 2)
    );
}

double BKllDecay::b(double q2) {
    return 2 * N(q2) * cache.m_l * beta_l(q2) * std::sqrt(lambda(q2)) * std::real(F_S(q2) * std::conj(F_V(q2)));
}

double BKllDecay::c(double q2) {
    return -0.25 * N(q2) * lambda(q2) * std::pow(beta_l(q2), 2) * (std::pow(std::abs(F_A(q2)), 2) + std::pow(std::abs(F_V(q2)), 2));
}

void BKllDecay::compute_binned_abc() {
    for (auto& v : cache.abc_binned) {
        v.clear();
    }
    cache.bin_widths.clear();

    if (!this->bins.has_value()) {
        LOG_WARN("BKllDecay::compute_binned_abc called without bins.");
        return;
    }

    constexpr double q2_min_eps = 1e-8;
    constexpr double endpoint_eps = 1e-2;
    constexpr double shat_eps = 1e-3;

    const double q2_qcdf_max = cache.m_b_PS * cache.m_b_PS * (1.0 - shat_eps);

    for (auto [q2_l, q2_u] : this->bins.value()) {
        const double low = std::max(q2_l, cache.q2_min + q2_min_eps);

        const double high = std::min({
            q2_u,
            cache.q2_max - endpoint_eps,
            q2_qcdf_max
        });

        if (q2_u > high) {
            // LOG_WARN(
            //     "BKll bin [", q2_l, ",", q2_u,
            //     "] clipped to [", low, ",", high,
            //     "] because q2_max = ", cache.q2_max,
            //     " and m_b_PS^2 safe max = ", q2_qcdf_max
            // );
        }

        if (!(low < high)) {
            LOG_WARN(
                "Skipping invalid BKll bin [",
                q2_l,
                ",",
                q2_u,
                "] after clipping to [",
                low,
                ",",
                high,
                "]"
            );

            cache.abc_binned[0].emplace_back(std::numeric_limits<double>::quiet_NaN());
            cache.abc_binned[1].emplace_back(std::numeric_limits<double>::quiet_NaN());
            cache.abc_binned[2].emplace_back(std::numeric_limits<double>::quiet_NaN());
            cache.bin_widths.emplace_back(std::numeric_limits<double>::quiet_NaN());
            continue;
        }

        const double width = high - low;

        try {
            cache.abc_binned[0].emplace_back(
                integrate([&](double q2) { return a(q2); }, low, high, 1e-3)
            );
            cache.abc_binned[1].emplace_back(
                integrate([&](double q2) { return b(q2); }, low, high, 1e-3)
            );
            cache.abc_binned[2].emplace_back(
                integrate([&](double q2) { return c(q2); }, low, high, 1e-3)
            );

            cache.bin_widths.emplace_back(width);
        } catch (const std::exception& e) {
            LOG_WARN(
                "BKll integration failed for bin [",
                q2_l,
                ",",
                q2_u,
                "] clipped to [",
                low,
                ",",
                high,
                "]: ",
                e.what()
            );

            cache.abc_binned[0].emplace_back(std::numeric_limits<double>::quiet_NaN());
            cache.abc_binned[1].emplace_back(std::numeric_limits<double>::quiet_NaN());
            cache.abc_binned[2].emplace_back(std::numeric_limits<double>::quiet_NaN());
            cache.bin_widths.emplace_back(std::numeric_limits<double>::quiet_NaN());
        }
    }
}

std::vector<ObservableValue> BKllDecay::dBR_dq2(Observables oid, bool br) {
    std::vector<ObservableValue> out;
    double br_factor = br ? cache.life_B : 1.0;

    for (size_t i = 0; i < this->bins.value().size(); i++) {
        const double integrated_rate = 2.0 * (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.0);
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

std::vector<ObservableValue> BKllDecay::A_FB(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = cache.abc_binned[1][i] / (2 * (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.)); 
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKllDecay::F_H(Observables oid) {
    std::vector<ObservableValue> out;
    for (size_t i = 0; i < this->bins.value().size(); i++) {
        double res = (cache.abc_binned[0][i] + cache.abc_binned[2][i]) / (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.); 
        out.emplace_back(ObservableMapper::to_id(oid), res, this->bins.value()[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKllDecay::Rm1_BK(Observables id, BKllConfig::B_Charge charge) {
    std::vector<ObservableValue> out;
    std::vector<double> Gamma_mu;
    std::vector<double> Gamma_e;

    set_lepton_gen_and_charge(BKllConfig::Lepton::MU, charge);

    for (size_t i = 0; i < this->bins.value().size(); i++)
        Gamma_mu.emplace_back(2 * (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.));

    set_lepton_gen_and_charge(BKllConfig::Lepton::E, charge);

    for (size_t i = 0; i < this->bins.value().size(); i++)
        Gamma_e.emplace_back(2 * (cache.abc_binned[0][i] + cache.abc_binned[2][i] / 3.));

    for (size_t i = 0; i < this->bins.value().size(); i++)
        out.emplace_back(ObservableMapper::to_id(id), Gamma_mu[i] / Gamma_e[i] - 1, this->bins.value()[i]);
    
    return out;
}

std::vector<ObservableValue> BKllDecay::compute_observable(Observables obs) {
    if (obs == Observables::R_1_B__K_L_L) {
        return Rm1_BK(obs, BKllConfig::B_Charge::B_PLUS);
    }

    if (obs == Observables::R_1_B0__K0_L_L) {
        return Rm1_BK(obs, BKllConfig::B_Charge::B_0);
    }

    auto it = BKllDecay::cfg_map.find(obs);
    if (it == BKllDecay::cfg_map.end()) {
        LOG_ERROR(
            "IndexError",
            "Observable",
            ObservableMapper::str(obs),
            "is not configured in BKllDecay::cfg_map"
        );
        return {};
    }

    auto flags = it->second;
    set_lepton_gen_and_charge(flags.first, flags.second);

    if (BKllDecay::dBR_dq2_ids.contains(ObservableMapper::to_id(obs))) return dBR_dq2(obs, true);
    if (BKllDecay::dG_dq2_ids.contains(ObservableMapper::to_id(obs)))  return dBR_dq2(obs, false);
    if (BKllDecay::A_FB_ids.contains(ObservableMapper::to_id(obs)))    return A_FB(obs);
    if (BKllDecay::F_H_ids.contains(ObservableMapper::to_id(obs)))     return F_H(obs);

    LOG_ERROR(
        "IndexError",
        "Observable",
        ObservableMapper::str(obs),
        "doesn't belong to the decay",
        DecayMapper::str(this->id)
    );

    return {};
}

std::vector<ObservableValue> BKllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}

void BKllDecay::set_n_threads(size_t n_threads) {
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
