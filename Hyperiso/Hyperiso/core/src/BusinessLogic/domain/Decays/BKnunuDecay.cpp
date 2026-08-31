#include "BKnunuDecay.h"
#include "wcoef_ids.hpp"

#include <algorithm>
#include <cmath>

namespace {
constexpr std::array<WCoef, 9> CNU_LEFT = {
    WCoef::CNU_L_EE, WCoef::CNU_L_EMU, WCoef::CNU_L_ETAU,
    WCoef::CNU_L_MUE, WCoef::CNU_L_MUMU, WCoef::CNU_L_MUTAU,
    WCoef::CNU_L_TAUE, WCoef::CNU_L_TAUMU, WCoef::CNU_L_TAUTAU
};
constexpr std::array<WCoef, 9> CNU_RIGHT = {
    WCoef::CNU_R_EE, WCoef::CNU_R_EMU, WCoef::CNU_R_ETAU,
    WCoef::CNU_R_MUE, WCoef::CNU_R_MUMU, WCoef::CNU_R_MUTAU,
    WCoef::CNU_R_TAUE, WCoef::CNU_R_TAUMU, WCoef::CNU_R_TAUTAU
};

// PDG 2022/23 value used in arXiv:2301.06990.  This becomes a FLIFE input in
// the follow-up assets patch so its uncertainty can participate in nuisances.
constexpr double TAU_TAU_SECONDS = 290.3e-15;

inline double kallen(double q2, double m1_sq, double m2_sq) {
    return q2 * q2 + m1_sq * m1_sq + m2_sq * m2_sq
         - 2.0 * (q2 * m1_sq + q2 * m2_sq + m1_sq * m2_sq);
}
}

void BKnunuDecay::load_params() {
    cache.G_F = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.lambda_t = (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE)
                   * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE));

    cache.m_Bp = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 521}, DataType::VALUE);
    cache.m_B0 = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 511}, DataType::VALUE);
    cache.m_Kp = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 321}, DataType::VALUE);
    cache.m_K0 = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 311}, DataType::VALUE);
    cache.tau_Bp = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 521}, DataType::VALUE) / HBAR;
    cache.tau_B0 = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 511}, DataType::VALUE) / HBAR;

    cache.m_tau = (*p)(ParamId{ParameterType::SM, "MASS", 15}, DataType::VALUE);
    cache.f_K = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {321, 1}}, DataType::VALUE);
    cache.f_B = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {521, 1}}, DataType::VALUE);
    cache.Vus = (*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE);
    cache.Vub = (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE);

    for (std::size_t i = 0; i < CNU_LEFT.size(); ++i) {
        cache.C_L[i] = w_proxy->getFR(WGroup::BNuNu, CNU_LEFT[i], w_config.order,
                                      ContributionType::TOTAL);
        cache.C_R[i] = w_proxy->getFR(WGroup::BNuNu, CNU_RIGHT[i], w_config.order,
                                      ContributionType::TOTAL);
    }
}

double BKnunuDecay::coefficient_sum_plus() const {
    double sum = 0.0;
    for (std::size_t i = 0; i < cache.C_L.size(); ++i) {
        sum += std::norm(cache.C_L[i] + cache.C_R[i]);
    }
    // Eq. (4) of arXiv:2301.06990 is already written after summing over the
    // three SM neutrino flavours and contains |C_L^SM|^2, not 3|C_L^SM|^2.
    // This 1/3 convention generalises it to flavour-diagonal and LFV WET input.
    return sum / 3.0;
}

double BKnunuDecay::f_plus(double q2, double m_B, double m_K) const {
    // Appendix A, Eqs. (A2), (A5), (A6), N=3, central values of Table VI.
    constexpr double a0 = 0.4742;
    constexpr double a1 = -0.894;
    constexpr double a2 = -0.44;
    constexpr double pole_mass = 5.4154;

    const double t_plus = std::pow(m_B + m_K, 2);
    const double t_zero = (m_B + m_K) * std::pow(std::sqrt(m_B) - std::sqrt(m_K), 2);
    const double z = (std::sqrt(t_plus - q2) - std::sqrt(t_plus - t_zero))
                   / (std::sqrt(t_plus - q2) + std::sqrt(t_plus - t_zero));

    // For N=3: sum_{n=0}^{2} a_n [z^n - (-1)^(n-3) n/3 z^3].
    const double z2 = z * z;
    const double z3 = z2 * z;
    const double series = a0 + a1 * (z - z3 / 3.0) + a2 * (z2 + 2.0 * z3 / 3.0);
    return series / (1.0 - q2 / (pole_mass * pole_mass));
}

double BKnunuDecay::loop_br(bool charged) {
    const double m_B = charged ? cache.m_Bp : cache.m_B0;
    const double m_K = charged ? cache.m_Kp : cache.m_K0;
    const double tau_B = charged ? cache.tau_Bp : cache.tau_B0;
    const double q2_max = std::pow(m_B - m_K, 2);
    const double csum = coefficient_sum_plus();

    const double pref = tau_B * std::pow(cache.G_F * cache.alpha_em, 2)
                      * std::norm(cache.lambda_t)
                      / (256.0 * std::pow(PI, 5) * std::pow(m_B, 3));

    auto differential = [&](double q2) {
        const double lambda = std::max(0.0, kallen(q2, m_B * m_B, m_K * m_K));
        const double fp = f_plus(q2, m_B, m_K);
        return pref * std::pow(lambda, 1.5) * csum * fp * fp;
    };

    double result = integrate(differential, 0.0, q2_max, 1e-6);
    // B0 -> K_S nu nu carries the |K_S> Clebsch-Gordan factor 1/sqrt(2).
    if (!charged) result *= 0.5;
    return result;
}

double BKnunuDecay::charged_tree_br() {
    // Eq. (20), narrow-width on-shell tau contribution.
    const double gamma_tau = HBAR / TAU_TAU_SECONDS;
    const double ckm2 = std::norm(cache.Vus * cache.Vub);
    const double mt2 = cache.m_tau * cache.m_tau;
    const double mK2 = cache.m_Kp * cache.m_Kp;
    const double mB2 = cache.m_Bp * cache.m_Bp;

    return cache.tau_Bp * std::pow(cache.G_F, 4) * ckm2
         * std::pow(cache.f_K * cache.f_B, 2)
         / (128.0 * PI2 * std::pow(cache.m_Bp, 3))
         * (cache.m_tau / gamma_tau)
         * std::pow(mt2 - mK2, 2)
         * std::pow(mB2 - mt2, 2);
}

std::vector<ObservableValue> BKnunuDecay::compute_observable(Observables obs) {
    double value = 0.0;
    switch (obs) {
        case Observables::BR_B__K_NU_NU:
            value = loop_br(true) + charged_tree_br();
            break;
        case Observables::BR_B0__KS_NU_NU:
            value = loop_br(false);
            break;
        default:
            LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs),
                      "doesn't belong to decay", DecayMapper::str(this->id));
    }
    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BKnunuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
