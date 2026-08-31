#include "BKstarnunuDecay.h"
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
constexpr double TAU_TAU_SECONDS = 290.3e-15;

inline double kallen(double q2, double m1_sq, double m2_sq) {
    return q2 * q2 + m1_sq * m1_sq + m2_sq * m2_sq
         - 2.0 * (q2 * m1_sq + q2 * m2_sq + m1_sq * m2_sq);
}
}

void BKstarnunuDecay::load_params() {
    cache.G_F = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.lambda_t = (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE)
                   * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE));

    cache.m_Bp = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 521}, DataType::VALUE);
    cache.m_B0 = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 511}, DataType::VALUE);
    cache.m_Kstarp = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 323}, DataType::VALUE);
    cache.m_Kstar0 = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", 313}, DataType::VALUE);
    cache.tau_Bp = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 521}, DataType::VALUE) / HBAR;
    cache.tau_B0 = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", 511}, DataType::VALUE) / HBAR;

    cache.m_tau = (*p)(ParamId{ParameterType::SM, "MASS", 15}, DataType::VALUE);
    cache.f_Kstar = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {323, 1}}, DataType::VALUE);
    cache.f_B = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {521, 1}}, DataType::VALUE);
    cache.Vus = (*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE);
    cache.Vub = (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE);

    cache.ff_charged = std::make_shared<BVFFCalculator>(521, 323, p, BV_FF_Src::BSZ_SR_LAT);
    cache.ff_neutral = std::make_shared<BVFFCalculator>(511, 313, p, BV_FF_Src::BSZ_SR_LAT);

    for (std::size_t i = 0; i < CNU_LEFT.size(); ++i) {
        cache.C_L[i] = w_proxy->getFR(WGroup::BNuNu, CNU_LEFT[i], w_config.order,
                                      ContributionType::TOTAL);
        cache.C_R[i] = w_proxy->getFR(WGroup::BNuNu, CNU_RIGHT[i], w_config.order,
                                      ContributionType::TOTAL);
    }
}

std::pair<double, double> BKstarnunuDecay::coefficient_sums() const {
    double minus = 0.0;
    double plus = 0.0;
    for (std::size_t i = 0; i < cache.C_L.size(); ++i) {
        minus += std::norm(cache.C_L[i] - cache.C_R[i]);
        plus += std::norm(cache.C_L[i] + cache.C_R[i]);
    }
    return {minus / 3.0, plus / 3.0};
}

double BKstarnunuDecay::loop_br(bool charged) {
    const double m_B = charged ? cache.m_Bp : cache.m_B0;
    const double m_V = charged ? cache.m_Kstarp : cache.m_Kstar0;
    const double tau_B = charged ? cache.tau_Bp : cache.tau_B0;
    BVFFCalculator& ff = *(charged ? cache.ff_charged : cache.ff_neutral);
    const auto [cminus, cplus] = coefficient_sums();

    const double mp = m_B + m_V;
    const double q2_max = std::pow(m_B - m_V, 2);
    // Eq. (7) of arXiv:2301.06990 contains (m_B + m_K*)^2 in the
    // numerator. F(q^2) in Eq. (8) is defined after factoring out this
    // quantity from the helicity amplitudes.
    const double pref = tau_B * std::pow(cache.G_F * cache.alpha_em, 2)
                      * std::norm(cache.lambda_t) * mp * mp
                      / (128.0 * std::pow(PI, 5) * std::pow(m_B, 3));

    // Eq. (6)-(8) of arXiv:2301.06990.  For C_R != 0 the axial A1/A12
    // helicities probe C_L-C_R while the vector V helicity probes C_L+C_R.
    // q^2 is multiplied through before evaluation, so the integrand remains
    // numerically finite at q^2=0 despite the conventional A12/q^2 term.
    auto differential = [&](double q2) {
        const double lambda = std::max(0.0, kallen(q2, m_B * m_B, m_V * m_V));
        const double sqrt_lambda = std::sqrt(lambda);
        const double A1 = ff.get(BV_FF::A1, q2);
        const double A12 = ff.get(BV_FF::A12, q2);
        const double V = ff.get(BV_FF::V, q2);

        const double axial = q2 * A1 * A1
                           + 32.0 * m_V * m_V * m_B * m_B / (mp * mp) * A12 * A12;
        const double vector = q2 * lambda / std::pow(mp, 4) * V * V;
        return pref * sqrt_lambda * (cminus * axial + cplus * vector);
    };

    return integrate(differential, 0.0, q2_max, 1e-6);
}

double BKstarnunuDecay::charged_tree_br() {
    // Eq. (22), narrow-width on-shell tau contribution.
    const double gamma_tau = HBAR / TAU_TAU_SECONDS;
    const double ckm2 = std::norm(cache.Vus * cache.Vub);
    const double mt2 = cache.m_tau * cache.m_tau;
    const double mV2 = cache.m_Kstarp * cache.m_Kstarp;
    const double mB2 = cache.m_Bp * cache.m_Bp;

    return cache.tau_Bp * std::pow(cache.G_F, 4) * ckm2
         * std::pow(cache.f_Kstar * cache.f_B, 2)
         / (128.0 * PI2 * std::pow(cache.m_Bp, 3))
         * (cache.m_tau / gamma_tau)
         * std::pow(mt2 - mV2, 2)
         * std::pow(mB2 - mt2, 2)
         * (1.0 + 2.0 * mV2 / mt2);
}

std::vector<ObservableValue> BKstarnunuDecay::compute_observable(Observables obs) {
    double value = 0.0;
    switch (obs) {
        case Observables::BR_B__KSTAR_NU_NU:
            value = loop_br(true) + charged_tree_br();
            break;
        case Observables::BR_B0__KSTAR0_NU_NU:
            value = loop_br(false);
            break;
        default:
            LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs),
                      "doesn't belong to decay", DecayMapper::str(this->id));
    }
    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BKstarnunuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
