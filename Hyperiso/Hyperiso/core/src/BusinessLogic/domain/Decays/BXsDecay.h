#ifndef BXSDECAY_H
#define BXSDECAY_H

#include "DecayParent.h"
#include "General.h"
#include <array>

/**
 * @brief Decay parent for the Bq > ll decays. Currently implements both CP-averaged and untagged Bs > mu+ mu- branching ratios and the CP-averaged Bd > mu+ mu- decays. 
 */
class BXsDecay : public DecayParent {

protected:
    static constexpr std::array<std::array<double, 8>, 8> gamma_0 {{
        {-4.0000, 2.6667,  0.0000, -0.2222, 0.0000,  0.0000, -0.8556,  1.0679},
        { 12.000, 0.0000,  0.0000,  1.3333, 0.0000,  0.0000,  5.1358,  2.5926},
        { 0.0000, 0.0000,  0.0000, -17.333, 0.0000,  2.0000, -2.1728,  0.5185},
        { 0.0000, 0.0000, -4.4444, -11.111, 0.4444,  0.8333, -0.6255, -3.6234},
        { 0.0000, 0.0000,  0.0000, -85.333, 0.0000,  20.000, -77.432,  244.30},
        { 0.0000, 0.0000, -28.444,  6.2222, 4.4444, -0.6667,  19.029,  58.914},
        { 0.0000, 0.0000,  0.0000,  0.0000, 0.0000,  0.0000,  10.667,  0.0000},
        { 0.0000, 0.0000,  0.0000,  0.0000, 0.0000,  0.0000, -3.5556,  9.3333}
    }};
    static constexpr std::array<double, 8> a_i {0.6086, 0.6956, 0.2609, -0.5217, 0.4086, -0.4230, -0.8994, 0.1456};
    static constexpr std::array<double, 8> d_i {1.4107, -0.8380, -0.4286, -0.0714, -0.6494, -0.0380, -0.0185, -0.0057};
    static constexpr std::array<std::pair<int, int>, 15> nonzero_K1 {{{0, 0}, {0, 1}, {0, 6}, {0, 7}, {1, 1}, {1, 6}, {1, 7}, {2, 6}, {3, 6}, {3, 7}, {4, 6}, {5, 6}, {6, 6}, {6, 7}, {7, 7}}};
    static constexpr std::array<std::pair<int, int>, 10> nonzero_K2 {{{0, 0}, {0, 1}, {0, 6}, {0, 7}, {1, 1}, {1, 6}, {1, 7}, {6, 6}, {6, 7}, {7, 7}}};

    template<std::size_t size>
    std::array<std::array<scalar_t, 8>, 8> fill_K(const std::vector<scalar_t>& flat_K, const std::array<std::pair<int, int>, size>& indices);

    double alpha_s_upsilon(double mb_1S);
    double delta(double E0, double mb_1S);
    double mc_muc(double mu_c);
    double mc_3GeV();
    double z(double mc_muc, double mb_1S);
    double Lb(double mu_b, double mb_1S);
    double Lc(double mu_c, double mc_muc);
    double LD(double Lb, double z);
    double a(double z);
    double b(double z);
    double r1(double az, double bz);
    double r2(double az, double bz);
    double r2_large_z(double z);
    double dr2_dlogz(double z);
    double r3(double a1, double b1);
    double r4(double bz, double a1, double b1);
    double r5(double a1, double b1);
    double r6(double az, double bz, double a1, double b1);
    double r8();
    complex_t G(double t);
    double phi_11(double phi_22); 
    double phi_12(double phi_22); 
    double phi_17(double phi_27); 
    double phi_18(double phi_27); 
    double phi_22(double delta, double z); 
    double phi_27(double delta, double z); 
    double phi_28(double phi_27); 
    double phi_47(double delta); 
    double phi_48(double phi_47); 
    double phi_77(double delta); 
    double phi_78(double delta); 
    double phi_88(double delta, double mb_1S, double ms); 
    double K_i7_1(int i, double ri, double Lb, double phi_i7);
    double K_77_1(double Lb, double phi_i7);
    double K_ij(int i, int j, double phi_ij);
    double r22(double z);
    double r22_large_z(double z);
    double h22(double z, double delta);
    double h27(double z, double delta);
    double h28(double z, double delta);
    double h88(double delta, double mb_1S, double ms);
    double F2nf(double z);
    double h77(double delta); 
    double phi_11_b0(double phi_22_b0);
    double phi_12_b0(double phi_22_b0);
    double phi_18_b0(double phi_28_b0);
    double phi_ij_b0(double phi_ij_1, double b0, double Lb, double hij);
    double K_17_b0(double K_27_b0);
    double K_27_b0(double b0, double r2, double az, double bz, double Lb, double phi_27_b0);
    double K_77_b0(double b0, double Lb, double phi_77_b0);
    double K_78_b0(double b0, double Lb);
    double F2a(double z);
    double F2na(double z);
    double phi_77_rem_int(double delta);
    double phi_77_rem(double phi_77_int, double delta, double alpha_upsilon);
    double K_11_rem(double K_22_rem);
    double K_12_rem(double K_22_rem);
    double K_17_rem(double K_27_rem, double K_78_1, double LD);
    double K_18_rem(double K_28_rem, double K_88_1, double LD);
    double K_27_large_z(double LD);
    double K_22_rem(double LD);
    double K_27_rem(double K_27_1, double K_77_1, double K_78_1, double K_47_1, double b0, double LD, double Lc, double Lb);
    double K_28_rem(double K_27_1, double K_78_1, double K_88_1, double K_48_1, double LD);
    double K_77_rem(double K_77_1, double phi_77_1, double phi_77_rem, double z, double LD, double Lb, double alpha_upsilon);
    double K_78_rem(double K_78_1, double LD);
    double K_88_rem(double K_88_1, double LD);
    double c(double z);
    double target(double p_22, double x1, double x2, double x5, double z);
    double x1();
    double x2();
    double x3(double P22_z0, double P22_z1, double c_z0, double c_z1);
    double x4(double P22_z0, double P22_z1, double c_z0, double c_z1);
    double x5(double P0, double K77_rem);
    double P0();
    double P11();
    double P12();

    template<std::size_t size>
    double gen_P00(const std::vector<scalar_t>& flat_K, const std::array<std::pair<int, int>, size>& indices);

    template<std::size_t size>
    double gen_P01(const std::vector<scalar_t>& flat_K, const std::array<std::pair<int, int>, size>& indices);

    double P22_rem(double x1, double x2, double x3, double x4, double x5, complex_t r21, double r22, double dr21_dlogz);
    double P(double alpha_mub, double p0, double p11, double p12, double p21, double p22, double p32);

    double g(double m_c_3gev, double m_b_1S);
    double C(double g, double m_b_mb, double m_c_3gev, double mu_G2, double rho_D3, double rho_LS3);

    double Kc(double eta);
    double Kt(double eta);
    double r(double m_b_mb, double m_b_1S);
    double N_eta_factor(double eta);
    double N(double Kc, double Kt, double r, double eta_factor, double lambda_2, double mc);

    double k_SL(double eta, double mu_W, double mu_b);
    double C2_em(double eta);
    double C8_em(double eta);
    complex_t C7_em(double eta, double C8_em, double C2_em);
    double epsilon_em(double inv_alpha_em, double alpha_mub, double C7_em, double k);

    double ckm(complex_t V_tb, complex_t V_ts, complex_t V_cb);
    double BR_B_Xs_gamma(double br_B__Xc_e_nu, double ckm, double inv_alpha_em, double C, double P, double N, double eps_em);

public:
    BXsDecay(QCDOrder order, double matching_scale, double hadronic_scale) {
        order = check_max_order(QCDOrder::NNLO);
        WilsonAdapter().build({WGroup::B, WGroup::BPrime}, matching_scale, hadronic_scale, order);
        this->w_proxy = std::make_shared<ObsWilsonProxy>();
        build_op_tree();
    }

    void build_op_tree() override;

};

#endif // __BXSDECAY_H__