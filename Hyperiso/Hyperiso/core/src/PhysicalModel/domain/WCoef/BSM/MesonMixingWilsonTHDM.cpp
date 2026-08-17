#include "MesonMixingWilsonTHDM.h"

#include <array>
#include <cmath>
#include <complex>
#include <unordered_set>
#include <stdexcept>

namespace {

struct THDMMixingResult {
    complex_t c1{0.0, 0.0};
    complex_t ct1{0.0, 0.0};
    complex_t c2{0.0, 0.0};
    complex_t ct2{0.0, 0.0};
    complex_t c3{0.0, 0.0};
    complex_t ct3{0.0, 0.0};
    complex_t c4{0.0, 0.0};
    complex_t c5{0.0, 0.0};
};

std::unordered_set<ParamId> thdm_mixing_sources(int light_down_pdg)
{
    return {
        {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},
        {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)},
        {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
        {ParameterType::SM, "MASS", 24},
        {ParameterType::BSM, "MASS", 37},
        {ParameterType::SM, "MASS", light_down_pdg},
        {ParameterType::SM, "MASS", 2},
        {ParameterType::SM, "GAUGE", 2},
        {ParameterType::BSM, "MINPAR", 3},
        {ParameterType::BSM, "MINPAR", 24},
        {ParameterType::SM, "VCKM", LhaID(0, 0)},
        {ParameterType::SM, "VCKM", LhaID(0, 1)},
        {ParameterType::SM, "VCKM", LhaID(0, 2)},
        {ParameterType::SM, "VCKM", LhaID(1, 0)},
        {ParameterType::SM, "VCKM", LhaID(1, 1)},
        {ParameterType::SM, "VCKM", LhaID(1, 2)},
        {ParameterType::SM, "VCKM", LhaID(2, 0)},
        {ParameterType::SM, "VCKM", LhaID(2, 1)},
        {ParameterType::SM, "VCKM", LhaID(2, 2)},
    };
}

THDMMixingResult charged_higgs_mixing(const ParamSrc& src,
                                      int light_down_pdg,
                                      int light_down_ckm_column)
{
    // Direct C++ translation of SuperIso's CM_calculator_chargedhiggs()
    // (core/Test/c_files_superiso/wilson2.c), with the original 1-based
    // arrays converted consistently to C++ 0-based arrays.
    const double m_h = src.get_val(ParameterType::BSM, "MASS", 37);
    const double m_w = src.get_val(ParameterType::SM, "MASS", 24);
    const double m_q = src.get_val(ParameterType::SM, "MASS", light_down_pdg);
    const double m_b = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1));
    const double g2 = src.get_val(ParameterType::SM, "GAUGE", 2);
    const double tan_beta = src.get_val(ParameterType::BSM, "MINPAR", 3);
    const int yukawa_type = static_cast<int>(std::lround(
        src.get_val(ParameterType::BSM, "MINPAR", 24)
    ));
    if (yukawa_type != 2) {
        throw std::runtime_error(
            "Native MesonMixingWilsonTHDM currently implements the SuperIso "
            "charged-Higgs formula for THDM type II only (MINPAR:24 = 2)."
        );
    }

    const double m_h2 = m_h * m_h;
    const double m_w2 = m_w * m_w;

    const std::array<double, 3> m_u{
        src.get_val(ParameterType::SM, "MASS", 2),
        src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 4),
        src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6),
    };
    const std::array<double, 3> m_u2{
        m_u[0] * m_u[0],
        m_u[1] * m_u[1],
        m_u[2] * m_u[2],
    };

    std::array<std::array<complex_t, 3>, 3> v{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            v[i][j] = src.get_val(ParameterType::SM, "VCKM", LhaID(i, j));
        }
    }

    THDMMixingResult result;
    const double inv_tb2 = std::pow(tan_beta, -2.0);
    const double inv_tb4 = std::pow(tan_beta, -4.0);
    const double tb2 = tan_beta * tan_beta;
    const double tb4 = tb2 * tb2;
    const double g24 = std::pow(g2, 4.0);
    const double mw4 = std::pow(m_w, 4.0);
    const double pi2 = std::pow(PI, 2.0);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const double d0h = D0(m_u2[i], m_u2[j], m_h2, m_w2);
            const double d0h_c = D0(m_u2[i], m_u2[j], m_h2, m_h2);
            const double d2h = D2p(m_u2[i], m_u2[j], m_h2, m_w2);
            const double d2h_c = D2p(m_u2[i], m_u2[j], m_h2, m_h2);

            // SuperIso 1-based: V[ie][3] V[je][3] conj(V[ie][gen]) conj(V[je][gen]).
            // Column 3 -> b -> 2 in 0-based indexing; gen=1 (Bd) / 2 (Bs)
            // -> columns 0 / 1 respectively.
            const complex_t ckm = v[i][2] * v[j][2]
                                * std::conj(v[i][light_down_ckm_column])
                                * std::conj(v[j][light_down_ckm_column]);

            result.c1 += g24 * ckm * m_u2[i] * m_u2[j]
                       * (2.0 * m_w2 * d0h * inv_tb2
                          - d2h_c * inv_tb4
                          - 2.0 * d2h * inv_tb2)
                       / (128.0 * pi2 * mw4);

            result.c2 += -g24 * (m_q * m_q) * ckm * m_u2[i] * m_u2[j]
                       * (d0h_c - 2.0 * d0h)
                       / (128.0 * pi2 * mw4);

            result.c4 += g24 * ckm
                       * (m_b * m_q * d2h * tb2 / m_w2
                          - m_b * m_q * m_u2[i] * m_u2[j]
                            * (d0h_c + d0h * (tb2 + inv_tb2))
                            / (4.0 * mw4))
                       / (16.0 * pi2);

            result.c5 += g24 * m_b * m_q * ckm * m_u2[i]
                       * (d2h_c - 2.0 * d2h)
                       / (32.0 * pi2 * mw4);

            result.ct1 += -g24 * (m_b * m_b) * (m_q * m_q) * ckm
                        * (d2h_c * tb4 + 2.0 * d2h * tb2)
                        / (128.0 * pi2 * mw4);

            result.ct2 += -g24 * (m_b * m_b) * ckm * m_u2[i] * m_u2[j]
                        * (d0h_c - 2.0 * d0h)
                        / (128.0 * pi2 * mw4);
        }
    }

    return result;
}

} // namespace

#define THDM_MATCHING_CTOR(CLASS, NAME, PDG, COMPUTE)                                  \
CLASS::CLASS() : WilsonCoefficient(NAME, GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) { \
    matching_info[QCDOrder::LO] = {                                                     \
        thdm_mixing_sources(PDG),                                                        \
        COMPUTE,                                                                         \
        get_lhaid_from_name(QCDOrder::LO)                                                \
    };                                                                                   \
}

// Bd: SuperIso gen=1 -> d column 0.
THDM_MATCHING_CTOR(C_mix_bd_1_THDM, "C_BD_1", 1, C_mix_bd_1_THDM::compute_LO)
complex_t C_mix_bd_1_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 1, 0).c1; }

THDM_MATCHING_CTOR(C_mix_bd_1_tilde_THDM, "CT_BD_1", 1, C_mix_bd_1_tilde_THDM::compute_LO)
complex_t C_mix_bd_1_tilde_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 1, 0).ct1; }

THDM_MATCHING_CTOR(C_mix_bd_2_THDM, "C_BD_2", 1, C_mix_bd_2_THDM::compute_LO)
complex_t C_mix_bd_2_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 1, 0).c2; }

THDM_MATCHING_CTOR(C_mix_bd_2_tilde_THDM, "CT_BD_2", 1, C_mix_bd_2_tilde_THDM::compute_LO)
complex_t C_mix_bd_2_tilde_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 1, 0).ct2; }

C_mix_bd_3_THDM::C_mix_bd_3_THDM()
    : WilsonCoefficient("C_BD_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_3_tilde_THDM::C_mix_bd_3_tilde_THDM()
    : WilsonCoefficient("CT_BD_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

THDM_MATCHING_CTOR(C_mix_bd_4_THDM, "C_BD_4", 1, C_mix_bd_4_THDM::compute_LO)
complex_t C_mix_bd_4_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 1, 0).c4; }

THDM_MATCHING_CTOR(C_mix_bd_5_THDM, "C_BD_5", 1, C_mix_bd_5_THDM::compute_LO)
complex_t C_mix_bd_5_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 1, 0).c5; }

// Bs: SuperIso gen=2 -> s column 1.
THDM_MATCHING_CTOR(C_mix_bs_1_THDM, "C_BS_1", 3, C_mix_bs_1_THDM::compute_LO)
complex_t C_mix_bs_1_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 3, 1).c1; }

THDM_MATCHING_CTOR(C_mix_bs_1_tilde_THDM, "CT_BS_1", 3, C_mix_bs_1_tilde_THDM::compute_LO)
complex_t C_mix_bs_1_tilde_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 3, 1).ct1; }

THDM_MATCHING_CTOR(C_mix_bs_2_THDM, "C_BS_2", 3, C_mix_bs_2_THDM::compute_LO)
complex_t C_mix_bs_2_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 3, 1).c2; }

THDM_MATCHING_CTOR(C_mix_bs_2_tilde_THDM, "CT_BS_2", 3, C_mix_bs_2_tilde_THDM::compute_LO)
complex_t C_mix_bs_2_tilde_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 3, 1).ct2; }

C_mix_bs_3_THDM::C_mix_bs_3_THDM()
    : WilsonCoefficient("C_BS_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_3_tilde_THDM::C_mix_bs_3_tilde_THDM()
    : WilsonCoefficient("CT_BS_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

THDM_MATCHING_CTOR(C_mix_bs_4_THDM, "C_BS_4", 3, C_mix_bs_4_THDM::compute_LO)
complex_t C_mix_bs_4_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 3, 1).c4; }

THDM_MATCHING_CTOR(C_mix_bs_5_THDM, "C_BS_5", 3, C_mix_bs_5_THDM::compute_LO)
complex_t C_mix_bs_5_THDM::compute_LO(const ParamSrc& src) { return charged_higgs_mixing(src, 3, 1).c5; }

#undef THDM_MATCHING_CTOR
