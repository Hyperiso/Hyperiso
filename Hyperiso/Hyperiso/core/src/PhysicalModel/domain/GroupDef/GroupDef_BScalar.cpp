#include "GroupDefinition.h"
#include "BScalarWilsonGroup.h"
#include "Math.h"
#include "ParameterProvider.h"

#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <string_view>

using CGS = CoefficientGroupSources;

namespace {
constexpr double kDefaultBsMassGeV = 5.36688;
constexpr double kDefaultAWidthGeV = 1.0e-6;

static double get_optional(const BlockSrc& src, std::string_view block, const LhaID& code, double fallback = 0.0) {
    try {
        return src.get_val(block, code);
    } catch (const std::exception&) {
        return fallback;
    }
}

static double get_optional(const BlockSrc& src, std::string_view block, int code, double fallback = 0.0) {
    return get_optional(src, block, LhaID(code), fallback);
}

static double get_optional(const BlockSrc& src, std::string_view block, std::initializer_list<int> code, double fallback = 0.0) {
    try {
        return src.get_val(block, code);
    } catch (const std::exception&) {
        return fallback;
    }
}

static bool almost_zero(double x) {
    return std::abs(x) < 1.0e-14;
}

static double chargino_mass(const BlockSrc& src, int i) {
    return src.get_val("WPARAM_SI_BSM", {13, i});
}

static double stop_mass(const BlockSrc& src, int a) {
    // WPARAM_SI_BSM {14,*} follows the helper ordering of the up-squark SLHA masses:
    // {1000002,1000004,1000006,2000002,2000004,2000006}.  The two stop entries
    // are therefore indices 2 and 5.
    return src.get_val("WPARAM_SI_BSM", {14, a == 0 ? 2 : 5});
}

static double stopmix(const BlockSrc& src, int mass_index, int gauge_index) {
    // SLHA STOPMIX is 1-based: STOPMIX(i,j), i,j = 1,2.
    return src.get_val("STOPMIX", {mass_index + 1, gauge_index + 1});
}

static double umix(const BlockSrc& src, int row, int col) {
    // SLHA UMIX is 1-based: UMIX(i,j), i,j = 1,2.
    return src.get_val("UMIX", {row + 1, col + 1});
}

static double vmix(const BlockSrc& src, int row, int col) {
    // SLHA VMIX is 1-based: VMIX(i,j), i,j = 1,2.
    return src.get_val("VMIX", {row + 1, col + 1});
}

static double nmamix(const BlockSrc& src, int a, int component) {
    // SLHA2 NMAMIX is 1-based: rows A_1,A_2 and columns Im(H_d), Im(H_u), Im(S).
    return src.get_val("NMAMIX", {a + 1, component + 1});
}

static double get_nmssm_optional(int nmssmrun_id, int extpar_id, double fallback = 0.0) {
    const ParameterProvider bsm(ParameterType::BSM);
    if (bsm.exists("NMSSMRUN", LhaID(nmssmrun_id))) {
        return static_cast<complex_t>(bsm("NMSSMRUN", LhaID(nmssmrun_id))).real();
    }

    const ParameterProvider passthrough(ParameterType::PASSTHROUGH);
    if (passthrough.exists("EXTPAR", LhaID(extpar_id))) {
        return static_cast<complex_t>(passthrough("EXTPAR", LhaID(extpar_id))).real();
    }

    return fallback;
}

static complex_t nmssm_pseudoscalar_threshold(const BlockSrc& src, double pseudoscalar_mass, int lepton_mass_slot) {
    const double lambda = get_nmssm_optional(1, 61, 0.0);
    const double kappa  = get_nmssm_optional(2, 62, 0.0);
    const double a_lambda = get_nmssm_optional(3, 63, 0.0);
    const double mu_eff = get_nmssm_optional(5, 65, get_optional(src, "HMIX", 1, 0.0));

    if (almost_zero(lambda) || almost_zero(mu_eff)) {
        return 0.0;
    }

    const double singlet_vev = mu_eff / lambda;
    const double denom = std::sqrt(2.0) * a_lambda + kappa * singlet_vev;
    if (almost_zero(denom)) {
        return 0.0;
    }

    const double tanb = src.get_val("HMIX", 2);
    const double mW = src.get_val("MASS", 24);
    const double mHpm = src.get_val("MASS", 37);
    const double mt_muW = src.get_val("WPARAM_MATCH_SM", 6);
    const double g2 = src.get_val("GAUGE", 2);
    const double gf = src.get_val("SMINPUTS", 2);
    const double v = std::sqrt(1.0 / (std::sqrt(2.0) * gf));

    const double v_delta_m_s = v / singlet_vev
        * (std::sqrt(2.0) * a_lambda - 2.0 * kappa * singlet_vev)
        / denom;

    complex_t cAH{0.0, -lambda * a_lambda / g2 / mW * tanb
                         * f30(mHpm * mHpm / mt_muW / mt_muW,
                               mW * mW / mt_muW / mt_muW)};

    complex_t cAc{0.0, 0.0};

    for (int i_chi = 0; i_chi < 2; ++i_chi) {
        for (int j_chi = 0; j_chi < 2; ++j_chi) {
            const double mi = chargino_mass(src, i_chi);
            const double mj = chargino_mass(src, j_chi);
            if (almost_zero(mi) || almost_zero(mj)) {
                continue;
            }

            const complex_t rA =
                -g2 / std::sqrt(2.0)
                    * (nmamix(src, 0, 0) * umix(src, 1, i_chi) * vmix(src, 0, j_chi)
                     + nmamix(src, 0, 1) * umix(src, 0, i_chi) * vmix(src, 1, j_chi))
                -lambda / std::sqrt(2.0)
                    * nmamix(src, 0, 2) * umix(src, 1, i_chi) * vmix(src, 1, j_chi);

            for (int a_stop = 0; a_stop < 2; ++a_stop) {
                const double mstop = stop_mass(src, a_stop);
                if (almost_zero(mstop)) {
                    continue;
                }

                const double tR = stopmix(src, a_stop, 1);
                const double tL = stopmix(src, a_stop, 0);
                const double g1_diag = (tR * tR - kron(a_stop + 1, 1))
                    * vmix(src, 1, j_chi) * umix(src, 1, i_chi)
                    - mt_muW / std::sqrt(2.0) / std::sin(std::atan(tanb)) / mW
                    * tL * tR * vmix(src, 1, j_chi) * umix(src, 1, i_chi);

                const double x_stop_i = std::pow(mstop / mi, 2.0);
                const double x_stop_j = std::pow(mstop / mj, 2.0);
                const double x_chi = std::pow(mi / mj, 2.0);

                const double loop_term =
                    v_delta_m_s * kron(i_chi, j_chi) * std::abs(mi / mW) * f80(x_stop_i)
                    - (rA * std::abs(mi / mj) * f30(x_stop_j, x_chi)
                       + rA * f40(x_stop_j, x_chi)).real();

                cAc = complex_t(cAc.real(), cAc.imag() + tanb / std::sqrt(2.0) * g1_diag * loop_term);
            }
        }
    }

    const complex_t cA = cAH + cAc;

    const double mBs = get_optional(src, "MASS", 531, kDefaultBsMassGeV);
    const double widthA = get_optional(src, "DECAY", LhaID(36), kDefaultAWidthGeV);
    const auto pole = complex_t{mBs * mBs - pseudoscalar_mass * pseudoscalar_mass,
                                pseudoscalar_mass * widthA};

    if (std::abs(pole) == 0.0) {
        return 0.0;
    }

    const double mb_muW_pole = src.get_val("WPARAM_MATCH_SM", {5, 2});
    const double sw2 = src.get_val("WPARAM_SI_SM", 4);
    const double ml_running = src.get_val("WPARAM_SI_SM", lepton_mass_slot);

    return (v_delta_m_s / 2.0) * mb_muW_pole / sw2 * ml_running * cA / pole;
}
} // namespace

static std::unordered_map<WCoefId, scalar_t>
BScalar_SUSY_Base1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src)
{
    auto out = BScalarCoefficientGroup::base_1_LO_calculation(coef_matching, src);

    const double mA1 = get_optional(src, "MASS", 36, 0.0);
    const double mA2 = get_optional(src, "MASS", 46, 0.0);
    const double h3  = get_optional(src, "MASS", 45, 0.0);

    // SLHA2/NMSSM: CP-even neutral Higgs masses are 25,35,45 and CP-odd
    // masses are 36,46.  If the NMSSM extension is absent, keep the generic
    // BScalar running result.
    if (almost_zero(h3) && almost_zero(mA2)) {
        return out;
    }

    const double mb_muW_pole = src.get_val("WPARAM_MATCH_SM", {5, 2});
    if (!almost_zero(mA1) && mA1 < mb_muW_pole) {
        for (int i = 0; i < 3; ++i) {
            out[WCoefMapper::to_id(WCoefMapper::cq2_for_lepton_index(i))] +=
                nmssm_pseudoscalar_threshold(src, mA1, WCoefMapper::lepton_mass_slot_from_index(i));
        }
    }

    return out;
}

static void Setup_BScalar_SUSY_Base1_LO(const BuildContext& ctx, CoefficientGroup& grp) {
    std::map<QCDOrder, CGS> m;

    CGS lo;
    lo.sources = {
        { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
                                   "WPARAM_RUN_SM", "WPARAM_SI_SM", "WPARAM_MATCH_SM", "WPARAM_SI_BSM" } },
        { ParameterType::SM,     { "SMINPUTS" } },
        { ParameterType::BSM,    { "GAUGE", "HMIX", "STOPMIX", "UMIX", "VMIX", "NMAMIX", "MASS" } }
    };

    // NMSSMRUN and EXTPAR are optional alternatives.  Do not inspect the
    // process-wide Parameters singleton from this build hook: group-definition
    // integration tests deliberately construct groups without a MemoryManager.
    // The low-mass NMSSM threshold reads whichever optional block is available
    // at evaluation time through ParameterProvider.

    lo.func = &BScalar_SUSY_Base1_LO_calculation;
    m[QCDOrder::LO] = lo;
    
    CGS nlo;
    nlo.sources = {
        { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
                                   "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
    };
    nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
    m[QCDOrder::NLO] = nlo;

    grp.add_sources(WilsonBasis::B_STANDARD, std::move(m));


}

namespace GroupDefinitions {
    const GroupDefinition& BScalar() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::BScalar);
            d.members = WCoefMapper::get_group(WGroup::BScalar);

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
            };
            lo.func = &BScalarCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            CGS nlo = lo; nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
            m[QCDOrder::NLO] = nlo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));

            d.setup[Model::SUSY].push_back(&Setup_BScalar_SUSY_Base1_LO);
            return d;
        }();
        return def;
    }
}
