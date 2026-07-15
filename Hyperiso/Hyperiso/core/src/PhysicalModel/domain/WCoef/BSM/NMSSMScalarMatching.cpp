#include "NMSSMScalarMatching.h"

#include <array>
#include <cmath>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>

#include "Logger.h"
#include "special_SUSY.h"
#include "scalar.h"
#include "ParameterProvider.h"

namespace nmssm_scalar_matching {
namespace {

constexpr double SQRT2 = 1.4142135623730950488;
constexpr double TINY = 1.0e-14;

using M22 = std::array<std::array<double, 2>, 2>;
using M23 = std::array<std::array<double, 3>, 2>;
using M33 = std::array<std::array<double, 3>, 3>;
using A222 = std::array<std::array<std::array<double, 2>, 2>, 2>;
using A322 = std::array<std::array<std::array<double, 2>, 2>, 3>;
using A333 = std::array<std::array<std::array<double, 3>, 3>, 3>;
using A3322 = std::array<std::array<std::array<std::array<double, 2>, 2>, 3>, 3>;

inline double real_value(const scalar_t& value) {
    return static_cast<complex_t>(value).real();
}

bool exists(ParameterType type, const std::string& block, const LhaID& id) {
    return ParameterProvider(type).exists(block, id);
}

double get(ParameterType type, const std::string& block, const LhaID& id) {
    return real_value(ParameterProvider(type)(block, id));
}

double get_required(ParameterType type, const std::string& block, const LhaID& id,
                    const std::string& label) {
    if (!exists(type, block, id)) {
        std::ostringstream os;
        os << "NMSSM scalar matching: missing " << label << " (" << block << ", " << id.to_string() << ")";
        throw std::runtime_error(os.str());
    }
    return get(type, block, id);
}

double get_nmssm_parameter(int nmssmrun_id, int extpar_id, const std::string& label,
                           bool allow_hmix_mu = false) {
    if (exists(ParameterType::BSM, "NMSSMRUN", LhaID(nmssmrun_id))) {
        return get(ParameterType::BSM, "NMSSMRUN", LhaID(nmssmrun_id));
    }
    if (exists(ParameterType::PASSTHROUGH, "EXTPAR", LhaID(extpar_id))) {
        return get(ParameterType::PASSTHROUGH, "EXTPAR", LhaID(extpar_id));
    }
    if (allow_hmix_mu && exists(ParameterType::BSM, "HMIX", LhaID(1))) {
        return get(ParameterType::BSM, "HMIX", LhaID(1));
    }
    std::ostringstream os;
    os << "NMSSM scalar matching: missing " << label
       << " (NMSSMRUN " << nmssmrun_id << " or EXTPAR " << extpar_id << ")";
    throw std::runtime_error(os.str());
}

inline double sqr(double x) { return x * x; }
inline double delta(int i, int j) { return i == j ? 1.0 : 0.0; }

void require_nonzero(double value, const std::string& label) {
    if (std::abs(value) < TINY) {
        throw std::runtime_error("NMSSM scalar matching: zero/near-zero " + label);
    }
}

} // namespace

bool is_active() {
    const ParameterProvider bsm(ParameterType::BSM);
    return (bsm.exists("MASS", LhaID(45)) && std::abs(real_value(bsm("MASS", LhaID(45)))) > TINY)
        || (bsm.exists("MASS", LhaID(46)) && std::abs(real_value(bsm("MASS", LhaID(46)))) > TINY);
}

Result compute(const ParamSrc& src, int lepton_mass_slot) {
    if (!is_active()) {
        return {};
    }

    const double lambda = get_nmssm_parameter(1, 61, "lambda");
    const double kappa_nmssm = get_nmssm_parameter(2, 62, "kappa");
    const double A_lambda = get_nmssm_parameter(3, 63, "A_lambda");
    const double lambda_s = get_nmssm_parameter(5, 65, "lambda*<S> = mu_eff", true);
    require_nonzero(lambda, "lambda");

    const double singlet_vev = lambda_s / lambda;
    require_nonzero(singlet_vev, "<S>");

    const double gf = get_required(ParameterType::SM, "SMINPUTS", LhaID(2), "G_F");
    const double v = std::sqrt(1.0 / std::sqrt(2.0) / gf);
    const double vdeltam_den = SQRT2 * A_lambda + kappa_nmssm * singlet_vev;
    require_nonzero(vdeltam_den, "sqrt(2) A_lambda + kappa <S>");
    const double v_deltam_s = v / singlet_vev
        * (SQRT2 * A_lambda - 2.0 * kappa_nmssm * singlet_vev) / vdeltam_den;

    const double mW = real_value(src.get_val(ParameterType::SM, "MASS", 24));
    const double mZ = get_required(ParameterType::SM, "SMINPUTS", LhaID(4), "M_Z");
    const double mHc = real_value(src.get_val(ParameterType::BSM, "MASS", 37));
    const double sw2 = real_value(src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", 4));
    const double mb_muW = real_value(src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}));
    const double mt_muW = real_value(src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", 6));
    const double g2 = real_value(src.get_val(ParameterType::SM, "GAUGE", 2));
    const double tanb = real_value(src.get_val(ParameterType::BSM, "MINPAR", 3));
    const double epsfac = real_value(src.get_val(ParameterType::WILSON, "EPSILON_SUSY", 5));
    const double ml = real_value(src.get_val(ParameterType::WILSON, "WPARAM_SI_SM", lepton_mass_slot));
    const double Qmatch = real_value(src.get_val(ParameterType::WILSON, "EW_SCALE", 1));

    require_nonzero(mW, "M_W");
    require_nonzero(mHc, "M_H+");
    require_nonzero(g2, "g2");
    require_nonzero(epsfac, "epsilon factor");

    const std::array<double, 3> mh = {
        get_required(ParameterType::SM, "MASS", LhaID(25), "MASS 25 (h1)"),
        get_required(ParameterType::BSM, "MASS", LhaID(35), "MASS 35 (h2)"),
        get_required(ParameterType::BSM, "MASS", LhaID(45), "MASS 45 (h3)")
    };
    const std::array<double, 2> ma = {
        get_required(ParameterType::BSM, "MASS", LhaID(36), "MASS 36 (a1)"),
        get_required(ParameterType::BSM, "MASS", LhaID(46), "MASS 46 (a2)")
    };
    const std::array<double, 3> mstop = {
        get_required(ParameterType::BSM, "MASS", LhaID(2000002), "MASS 2000002 (u_R-like squark)"),
        get_required(ParameterType::BSM, "MASS", LhaID(1000006), "MASS 1000006 (stop 1)"),
        get_required(ParameterType::BSM, "MASS", LhaID(2000006), "MASS 2000006 (stop 2)")
    };
    const double m_snutau = get_required(ParameterType::BSM, "MASS", LhaID(1000016), "MASS 1000016 (tau sneutrino)");
    const double Au = get_required(ParameterType::BSM, "AU", LhaID(1, 1), "AU(1,1)");

    std::array<double, 2> mch {};
    for (int j = 0; j < 2; ++j) {
        mch[j] = real_value(src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", {13, j}));
        require_nonzero(mch[j], "chargino mass");
    }

    M22 U {}, V {}, stopmix {};
    M33 H {}, TU {};
    M23 A {};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            U[i][j] = get_required(ParameterType::BSM, "UMIX", LhaID(i + 1, j + 1), "UMIX");
            V[i][j] = get_required(ParameterType::BSM, "VMIX", LhaID(i + 1, j + 1), "VMIX");
            stopmix[i][j] = get_required(ParameterType::BSM, "STOPMIX", LhaID(i + 1, j + 1), "STOPMIX");
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int c = 0; c < 3; ++c) {
            H[a][c] = get_required(ParameterType::BSM, "NMHMIX", LhaID(a + 1, c + 1), "NMHMIX");
        }
    }
    for (int a = 0; a < 2; ++a) {
        for (int c = 0; c < 3; ++c) {
            A[a][c] = get_required(ParameterType::BSM, "NMAMIX", LhaID(a + 1, c + 1), "NMAMIX");
        }
    }

    TU[0][0] = 1.0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            TU[i + 1][j + 1] = stopmix[i][j];
        }
    }

    const double beta = std::atan(tanb);
    const double sinb = std::sin(beta);
    const double cosb = std::cos(beta);
    require_nonzero(sinb, "sin(beta)");
    require_nonzero(cosb, "cos(beta)");
    const double vu = std::sqrt(sqr(sinb) / std::sqrt(2.0) / gf);
    const double vd = vu / tanb;

    A222 R {};
    A322 Q {};
    A3322 G1 {};
    A333 T2 {};
    std::array<std::array<std::array<double, 3>, 3>, 2> T1 {};

    for (int a = 0; a < 2; ++a) {
        for (int l = 0; l < 2; ++l) {
            for (int j = 0; j < 2; ++j) {
                R[a][l][j] = -g2 / SQRT2
                    * (A[a][0] * U[1][l] * V[1][j] + A[a][1] * U[0][l] * V[1][j])
                    - lambda / SQRT2 * A[a][2] * U[1][l] * V[1][j];
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int l = 0; l < 2; ++l) {
            for (int j = 0; j < 2; ++j) {
                Q[a][l][j] = g2 / SQRT2
                    * (H[a][0] * U[1][l] * V[1][j] + H[a][1] * U[0][l] * V[1][j])
                    - lambda / SQRT2 * H[a][2] * U[1][l] * V[1][j];
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 2; ++j) {
                for (int l = 0; l < 2; ++l) {
                    G1[i][k][j][l] = (TU[i][1] * TU[k][1] - delta(i, 0) * delta(k, 0)) * V[0][l] * U[1][j]
                        - mt_muW / SQRT2 / sinb / mW * TU[i][2] * TU[k][1] * V[1][l] * U[1][j];
                }
            }
        }
    }
    for (int a = 0; a < 2; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                T1[a][i][k] = (TU[i][2] * TU[k][1] - TU[i][1] * TU[k][2])
                    * ((lambda / SQRT2 * (vd * A[a][2] + singlet_vev * A[a][0])) - Au * A[a][1]);
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                T2[a][i][k] = -mt_muW / (2.0 * mW)
                    * (2.0 * mt_muW * H[a][1] * (TU[i][1] * TU[k][1] + TU[i][2] * TU[k][2])
                       + ((lambda / SQRT2 * (vd * H[a][2] + singlet_vev * H[a][0])) + Au * H[a][1])
                           * (TU[i][2] * TU[k][1] + TU[i][1] * TU[k][2]))
                    + mZ / 2.0 / std::sqrt(1.0 - sw2) * (1.0 - 4.0 / 3.0 * sw2) * H[a][1]
                        * (TU[i][0] * TU[k][0] + TU[i][1] * TU[k][1])
                    + 2.0 / 3.0 * mW * sw2 / (1.0 - sw2) * H[a][1] * TU[i][2] * TU[k][2];
            }
        }
    }

    complex_t cq1_h = 0.0;
    for (int a = 0; a < 3; ++a) {
        require_nonzero(mh[a], "CP-even Higgs mass");
        cq1_h += (sqr(mHc / mW) * sqr(H[a][0]) * f30(sqr(mHc / mt_muW), sqr(mW / mt_muW))
                  + sqr(mt_muW) * sqr(mh[a]) / sqr(mW) / sqr(mHc)
                      * f30(sqr(mt_muW / mHc), sqr(mt_muW / mW))) / sqr(mh[a]);
    }
    cq1_h *= -ml / 4.0 * sqr(tanb);

    complex_t cq2_h = 0.0;
    for (int a = 0; a < 2; ++a) {
        require_nonzero(ma[a], "CP-odd Higgs mass");
        cq2_h += ((sqr(mHc / mW) * sqr(A[a][0]) + delta(a, 1) * A[a][0])
                      * f30(sqr(mHc / mt_muW), sqr(mW / mt_muW))
                  + sqr(mt_muW) * sqr(ma[a]) / sqr(mW) / sqr(mHc)
                      * f30(sqr(mt_muW / mHc), sqr(mt_muW / mW))) / sqr(ma[a]);
    }
    cq2_h *= ml / 4.0 * sqr(tanb);

    const complex_t I(0.0, 1.0);
    const complex_t ca_h = -I * lambda * A_lambda / g2 / mW * tanb
        * f30(sqr(mHc / mt_muW), sqr(mW / mt_muW));

    complex_t cq1_c = 0.0;
    for (int a = 0; a < 3; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                for (int j = 0; j < 2; ++j) {
                    for (int l = 0; l < 2; ++l) {
                        cq1_c += G1[i][k][j][l] / sqr(mh[a]) * (
                            SQRT2 * sqr(H[a][0]) * mch[j] / mW / cosb * delta(i, k) * delta(l, j)
                                * f80(sqr(mstop[i] / mch[j]))
                            - 2.0 * SQRT2 * H[a][0] / g2 * delta(i, k)
                                * (Q[a][l][j] * f40(sqr(mstop[i] / mch[l]), sqr(mch[j] / mch[l]))
                                   + mch[j] / mch[l] * Q[a][j][l]
                                       * f30(sqr(mstop[i] / mch[l]), sqr(mch[j] / mch[l])))
                            + 2.0 * SQRT2 * H[a][0] * T2[a][i][k] * mch[j] / sqr(mstop[k]) * delta(l, j)
                                * f30(sqr(mstop[i] / mstop[k]), sqr(mch[j] / mstop[k]))
                            + sqr(mh[a] / mch[j]) * delta(i, k)
                                * (U[1][j] * V[0][l]
                                       * f50(sqr(mstop[i] / mch[j]), sqr(mch[l] / mch[j]), sqr(m_snutau / mch[l]))
                                   - mch[l] / mch[j] * U[1][l] * V[0][j]
                                       * f60(sqr(mstop[i] / mch[j]), sqr(mch[l] / mch[j]), sqr(m_snutau / mch[l])))
                        );
                    }
                }
            }
        }
    }
    cq1_c *= ml / 4.0 * sqr(tanb);

    complex_t cq2_c = 0.0;
    for (int a = 0; a < 2; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                for (int j = 0; j < 2; ++j) {
                    for (int l = 0; l < 2; ++l) {
                        cq2_c += G1[i][k][j][l] / sqr(ma[a]) * (
                            SQRT2 * sqr(A[a][0]) * mch[j] / mW / cosb * delta(i, k) * delta(l, j)
                                * f80(sqr(mstop[i] / mch[j]))
                            - 2.0 * SQRT2 * A[a][0] / g2 * delta(i, k)
                                * (-R[a][l][j] * f40(sqr(mstop[i] / mch[l]), sqr(mch[j] / mch[l]))
                                   + mch[j] / mch[l] * R[a][j][l]
                                       * f30(sqr(mstop[i] / mch[l]), sqr(mch[j] / mch[l])))
                            - SQRT2 * A[a][0] * T1[a][i][k] * mt_muW * mch[j] / sqr(mstop[k]) * delta(l, j)
                                * f30(sqr(mstop[i] / mstop[k]), sqr(mch[j] / mstop[k]))
                            + sqr(ma[a] / mch[j]) * delta(i, k)
                                * (U[1][j] * V[0][l]
                                       * f50(sqr(mstop[i] / mch[j]), sqr(mch[l] / mch[j]), sqr(m_snutau / mch[l]))
                                   - mch[l] / mch[j] * U[1][l] * V[0][j]
                                       * f60(sqr(mstop[i] / mch[j]), sqr(mch[l] / mch[j]), sqr(m_snutau / mch[l])))
                        );
                    }
                }
            }
        }
    }
    cq2_c *= -ml / 4.0 * sqr(tanb);

    complex_t ca_c = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int l = 0; l < 2; ++l) {
                ca_c += I * tanb / SQRT2 * G1[i][i][j][l]
                    * (v_deltam_s * delta(l, j) * std::abs(mch[j] / mW) * f80(sqr(mstop[i] / mch[j]))
                       - (R[0][j][l] * std::abs(mch[j] / mch[l])
                              * f30(sqr(mstop[i] / mch[l]), sqr(mch[j] / mch[l]))
                          - R[0][l][j] * f40(sqr(mstop[i] / mch[l]), sqr(mch[j] / mch[l]))));
            }
        }
    }

    complex_t cq1 = (cq1_h + cq1_c) * mb_muW / sw2 / epsfac;
    complex_t cq2 = (cq2_h + cq2_c) * mb_muW / sw2 / epsfac;
    const complex_t ca = ca_h + ca_c;

    // This is the first SuperIso threshold term.  The lower-scale cases are
    // handled in GroupDef_BScalar so the contribution is not counted twice.
    if (ma[0] > Qmatch) {
        cq2 += -v_deltam_s / 2.0 * mb_muW / sw2 * ml * ca / sqr(ma[0]);
    }

    LOG_INFO("NMSSM scalar matching enabled:",
             "lambda=", lambda, "kappa=", kappa_nmssm,
             "A_lambda=", A_lambda, "mu_eff=", lambda_s,
             "M_h3=", mh[2], "M_a2=", ma[1]);

    return {true, scalar_t(cq1), scalar_t(cq2)};
}

} // namespace nmssm_scalar_matching
