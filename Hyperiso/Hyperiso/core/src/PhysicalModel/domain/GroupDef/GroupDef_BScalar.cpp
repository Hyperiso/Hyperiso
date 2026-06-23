#include "GroupDefinition.h"
#include "BScalarWilsonGroup.h"
#include "Math.h"

#include <array>
#include <cmath>
#include <exception>
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

static complex_t nmssm_pseudoscalar_threshold(const BlockSrc& src, double pseudoscalar_mass) {
    const double lambda = get_optional(src, "EXTPAR", 61, 0.0);
    const double kappa  = get_optional(src, "EXTPAR", 62, 0.0);
    const double a_lambda = get_optional(src, "EXTPAR", 63, 0.0);
    const double mu_eff = get_optional(src, "EXTPAR", 65, get_optional(src, "HMIX", 1, 0.0));

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
    const double mb_running = src.get_val("WPARAM_SI_SM", 3);

    return (v_delta_m_s / 2.0) * mb_muW_pole / sw2 * mb_running * cA / pole;
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
        out[WCoefMapper::to_id(WCoef::CQ2)] += nmssm_pseudoscalar_threshold(src, mA1);
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
        { ParameterType::BSM,    { "GAUGE", "HMIX", "STOPMIX", "UMIX", "VMIX", "NMAMIX", "MASS", "EXTPAR" } }
    };
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


// #include "GroupDefinition.h"
// #include "BScalarWilsonGroup.h"
// #include "Math.h"

// using CGS = CoefficientGroupSources;

// static std::unordered_map<WCoefId, scalar_t>
// BScalar_SUSY_Base1_LO_calculation(
//     const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
//     const BlockSrc& src)
// {
//     const auto ids = WCoefMapper::get_group(WGroup::BScalar);

//     const auto itLO = coef_matching.find(QCDOrder::LO);
//     const auto* matchLO = (itLO != coef_matching.end()) ? &itLO->second : nullptr;

//     auto getM = [&](WCoefId c) -> scalar_t {
//         if (!matchLO) return scalar_t(0);
//         auto it = matchLO->find(c);
//         return (it != matchLO->end()) ? it->second : scalar_t(0);
//     };

//     const double eta    = src.get_val("WPARAM_RUN_SM",2);
//     const double beta_0 = src.get_val("WPARAM_SI_SM",5);

//     const double fact = std::pow(eta, -4.0 / beta_0);

//     std::unordered_map<WCoefId, scalar_t> out;
//     out.reserve(ids.size());
//     for (auto c : ids) {
//         out[WCoefMapper::to_id(c)] = fact * getM(WCoefMapper::to_id(c));
//     }

//     complex_t coeff_temp2 = 0;
//     if(src.get_val("MASS",46)!=0.||src.get_val("MASS",45)!=0.) {
//             double mass_b_2 = src.get_val("QCD",LhaID(5, 2));
//             double mA = src.get_val("MASS",36);
// 			if(mA < mass_b_2) {	
//             double lambdaNMSSM = 1;
//             double lambdaSNMSSM = 1;
//             double AlambdaNSSM = 1;
//             double kappaNMSSM = 1;
//             double m_Bs = 1;
//             // double mass_nutl = 1;
            
//             double sw2 = src.get_val("WPARAM_SI_SM",4);
//             double mH = src.get_val("MASS",37);
		    
//             double mass_top_muW = src.get_val("WPARAM_MATCH_SM",6);
//             double g2 = src.get_val("GAUGE",2);
//             double tanb = src.get_val("HMIX",2);
//             double mW = src.get_val("MASS",24);

//             // double mH0[4],mA0[3],mstop[3];
//             double mstop[3];

//             mstop[0]=src.get_val("MASS",2000002); //mass upr, is that right ?
//             mstop[1]=src.get_val("MASS",1000006);
//             mstop[2]=src.get_val("MASS",2000006);

//             complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/g2/mW*tanb*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)};
//             complex_t CAc{};
//             double s=lambdaSNMSSM/lambdaNMSSM;
//             double v=sqrt(1./sqrt(2.)/src.get_val("SMINPUTS",2));
//             double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);

//             // double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
//             double Ralj[3][3][3],G1[4][4][3][3];
//             // double T2[4][4][4];
//             std::array<std::array<double,4>,4> TU;
//             // double vu=sqrt(pow(sin(atan(tanb)),2.)/sqrt(2.)/src.get_val("SMINPUTS",2));
//             // double vd=vu/tanb;

//             TU[1][1]=1.;
//             for(int ie=0;ie<2;ie++){
//                 for(int je=0;je<2;je++) {
//                     TU[ie+1][je+1]=src.get_val("STOPMIX",{ie+1, je+1});
//                 }
//             }

//             for(int je=0;je<2;je++) {
//                 for(int le=0;le<2;le++) {
//                     for(int ae=0;ae<3;ae++) {
//                         if (ae <3 ){
//                             Ralj[ae][le][je]=-g2/sqrt(2.)*(src.get_val("NMAMIX",{ae+1, 1+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1})+src.get_val("NMAMIX",{ae+1, 2+1})*src.get_val("UMIX",{1+1, le+1})*src.get_val("VMIX",{2+1, je+1}))-lambdaNMSSM/sqrt(2.)*src.get_val("NMAMIX",{ae+1, 3+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1});
//                         }
//                         // Qalj[ae][le][je]=g2/sqrt(2.)*(src.get_val("NMHMIX",{ae+1,1+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1})+src.get_val("NMHMIX",{ae+1, 2+1})*src.get_val("UMIX",{1+1, le+1})*src.get_val("VMIX",{2+1, je+1}))-lambdaNMSSM/sqrt(2.)*src.get_val("NMHMIX",{ae+1, 3+1})*src.get_val("UMIX",{2+1, le+1})*src.get_val("VMIX",{2+1, je+1});
//                         for(int ke=1;ke<=3;ke++) {
//                             G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*src.get_val("VMIX",{1+1, le+1})*src.get_val("UMIX",{2+1, je+1})-mass_top_muW/sqrt(2.)/sin(atan(tanb))/mW*TU[ae][3]*TU[ke][2]*src.get_val("VMIX",{2+1, le+1})*src.get_val("UMIX",{2+1, je+1});
//                         }
//                     }
//                 }
//             }
//             for(int ae=0;ae<3;ae++) {
//                 for(int je=0;je<2;je++) {
//                     for(int le=0;le<2;le++) {
//                         CAc = complex_t(CAc.real(), CAc.imag()+(tanb)/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs(src.get_val("WPARAM_SI_BSM",je)/mW)*f80(pow(mstop[ae-1]/src.get_val("WPARAM_SI_BSM",je),2.))-(Ralj[1][je][le]*fabs(src.get_val("WPARAM_SI_BSM",je)/src.get_val("WPARAM_SI_BSM",le))*f30(pow(mstop[ae-1]/src.get_val("WPARAM_SI_BSM",le),2.),pow(src.get_val("WPARAM_SI_BSM",je)/src.get_val("WPARAM_SI_BSM",le),2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/src.get_val("WPARAM_SI_BSM",le),2.),pow(src.get_val("WPARAM_SI_BSM",je)/src.get_val("WPARAM_SI_BSM",le),2.)))));
//                     }
//                 }
//             }
//             complex_t CA=CAH+CAc;
//             double width_A0=1.e-6;
//             auto denom = complex_t{ m_Bs*m_Bs - mA*mA, mA*width_A0 };

//             auto num = v_deltam_s/2. * mass_b_2 / sw2
//                     * src.get_val("WPARAM_SI_SM", 3)
//                     * CA;

//             coeff_temp2 += num / denom;
//             // coeff_temp2+=complex_t{v_deltam_s/2.*mass_b_2/sw2*src.get_val("WPARAM_SI_SM",3)*CA/(m_Bs*m_Bs-mA*mA,mA*width_A0)};
//         }
//     }
//     out[WCoefMapper::to_id(WCoef::CQ2)] += coeff_temp2;
//     return out;
// }

// static void Setup_BScalar_SUSY_Base1_LO(const BuildContext& ctx, CoefficientGroup& grp) {
//     std::map<QCDOrder, CGS> m;

//     CGS lo;
//     lo.sources = {
//         { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
//                                    "WPARAM_RUN_SM", "WPARAM_SI_SM", "WPARAM_MATCH_SM" } },
//         { ParameterType::SM,     { "MASS", "SMINPUTS" } },
//         { ParameterType::BSM,    { "GAUGE", "HMIX", "STOPMIX", "UMIX", "VMIX", "NMAMIX", "NMHMIX", "MASS" } }
//     };
//     lo.func = &BScalar_SUSY_Base1_LO_calculation;
//     m[QCDOrder::LO] = lo;
    
//     CGS nlo;
//     nlo.sources = {
//         { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
//                                    "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
//     };
//     nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
//     m[QCDOrder::NLO] = nlo;

//     grp.add_sources(WilsonBasis::B_STANDARD, std::move(m));


// }

// namespace GroupDefinitions {
//     const GroupDefinition& BScalar() {
//         static const GroupDefinition def = []{
//             GroupDefinition d;
//             d.id = GroupMapper::to_id(WGroup::BScalar);
//             d.members = WCoefMapper::get_group(WGroup::BScalar);

//             std::map<QCDOrder, CGS> m;
//             CGS lo;
//             lo.sources = {
//                 { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
//             };
//             lo.func = &BScalarCoefficientGroup::base_1_LO_calculation;
//             m[QCDOrder::LO] = lo;

//             CGS nlo = lo; nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
//             m[QCDOrder::NLO] = nlo;

//             d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));

//             d.setup[Model::SUSY].push_back(&Setup_BScalar_SUSY_Base1_LO);
//             return d;
//         }();
//         return def;
//     }
// }
