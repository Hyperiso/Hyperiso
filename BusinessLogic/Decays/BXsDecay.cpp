#include "BXsDecay.h"

template<std::size_t size>
std::array<std::array<scalar_t, 8>, 8> BXsDecay::fill_K(const std::vector<scalar_t>& flat_K, const std::array<std::pair<int, int>, size>& indices) {
    std::array<std::array<scalar_t, 8>, 8> K {};
    int k {0};
    for (auto [i, j] : indices) {
        K[j][i] = flat_K[k];
        K[i][j] = flat_K[k++];
    }
    return K;
}

double BXsDecay::alpha_s_mub() {
    return QCDHelper::alpha_s(winfo.hadronic_scale);
}

double BXsDecay::alpha_s_upsilon(double mb_1S) {
    return QCDHelper::alpha_s(mb_1S);
}

double BXsDecay::eta(double alpha_mub) {
    return QCDHelper::alpha_s(winfo.matching_scale) / alpha_mub;
}

double BXsDecay::delta(double E0, double mb_1S) {
    return 1 - 2 * E0 / mb_1S;
}

double BXsDecay::mc_muc(double mu_c) {
    return QCDHelper::msbar_mass(4, mu_c, "running");;
}

double BXsDecay::mc_3GeV() {
    return QCDHelper::msbar_mass(4, 3, "running");
}

double BXsDecay::z(double mc_muc, double mb_1S) {
    return std::pow(mc_muc / mb_1S, 2);
}

double BXsDecay::Lb(double mb_1S) {
    return 2 * std::log(this->winfo.hadronic_scale / mb_1S);
}

double BXsDecay::Lc(double mu_c, double mc_muc) {
    return 2 * std::log(mu_c / mc_muc);
}

double BXsDecay::LD(double Lb, double z) {
    return Lb - std::log(z);
}

double BXsDecay::a(double z) {
    if(fpeq(z, 1.)) return 4.0859;
	if(fpeq(z, 0.)) return 0.;
    if(std::abs(z) > 0.4) LOG_WARN("The value of z in BXsDecay::a(double z) shouldn't exceed 0.4.");
		
	double Lz = std::log(z);
    double Lz2 = Lz * Lz;

	return 16./9.*((5./2.-PI2/3.-3.*ZETA3+(5./2.-3./4.*PI2)*Lz+Lz2/4.+Lz2*Lz/12.)*z
	+(7./4.+2./3.*PI2-0.5*PI2*Lz-Lz2/4.+Lz2*Lz/12.)*std::pow(z,2)+(-7./6.-PI2/4.+2.*Lz-3./4.*Lz2)*std::pow(z,3)
	+(457./216.-5./18.*PI2-Lz/72.-5./6.*Lz2)*std::pow(z,4)+(35101./8640.-35./72.*PI2-185./144.*Lz-35./24.*Lz2)*std::pow(z,5)
	+(67801./8000.-21./20.*PI2-3303./800.*Lz-63./20.*Lz2)*std::pow(z,6));
}

double BXsDecay::b(double z) {
    if(fpeq(z, 1.)) return 0.0316;
	if(fpeq(z, 0.)) return 0.;
    if(std::abs(z) > 0.4) LOG_WARN("The value of z in BXsDecay::b(double z) shouldn't exceed 0.4.");
		
	double Lz = std::log(z);
    double Lz2 = Lz * Lz;

	return -8./9.*((-3.+PI2/6.-Lz)*z-2./3.*PI2*std::pow(z,1.5)+(0.5+PI2-2.*Lz-0.5*Lz2)*z*z
	+(-25./12.-PI2/9.-19./18.*Lz+2.*Lz2)*z*z*z+(-1376./225.+137./30.*Lz+2.*Lz2+2./3.*PI2)*std::pow(z,4)
	+(-131317./11760.+887./84.*Lz+5.*Lz2+5./3.*PI2)*std::pow(z, 5)+(-2807617./97200.+16597./540.*Lz+14.*Lz2+14./3.*PI2)*std::pow(z,6));
}

double BXsDecay::r1(double az, double bz) {
    return 833. / 729 - (az + bz) / 3;
}

double BXsDecay::r2(double az, double bz) {
    return -1666. / 243 + 2 * (az + bz);
}

double BXsDecay::r2_large_z(double z) {
    return -1666. / 243 + 2 * (104 * std::log(z) + 314) / 81;
}

double BXsDecay::dr2_dlogz(double z) {
    double lz = std::log(z);
    double lz2 = lz * lz;
    double lz3 = lz2 * lz;
    return                   z * (224. / 9           -  112 * PI2 / 27 - 32 * ZETA3 / 3 + (112. / 9 - 8 * PI2 / 3) * lz +   16 * lz2 / 9 +  8 * lz3 / 27 )
            +   std::pow(z, 2) * (128. / 9           -   16 * PI2 / 27                  + (64. / 9 - 32 * PI2 / 9) * lz +    8 * lz2 / 9 + 16 * lz3 / 27 )
            + std::pow(z, 1.5) * (                       16 * PI2 / 9                                                                                    )
            +   std::pow(z, 3) * (620. / 81          -   56 * PI2 / 27                  + 392. / 27                * lz -   56 * lz2 / 3                 )
            +   std::pow(z, 4) * (397372. / 6075     -  704 * PI2 / 81                  - 18512. / 405             * lz -  704 * lz2 / 27                )
            +   std::pow(z, 5) * (-1199585. / 23814  - 1900 * PI2 / 81                  - 82130. / 567             * lz - 1900 * lz2 / 27                )
            +   std::pow(z, 6) * (-17917342. / 91125 - 3248 * PI2 / 45                  - 988402. / 2025           * lz - 3248 * lz2 / 15                );
}

double BXsDecay::r3(double a1, double b1) {
    double Xb = -0.1684;
    return 2392. / 243 + 8 * PI / (3 * std::sqrt(3)) + 32 * Xb / 9 - a1 + 2 * b1;
}

double BXsDecay::r4(double bz, double a1, double b1) {
    double Xb = -0.1684;
    return -761. / 729 - 4 * PI / (9 * std::sqrt(3)) - 16 * Xb / 27 + a1 / 6 + 5 * b1 /3 + 2 * bz;
}

double BXsDecay::r5(double a1, double b1) {
    double Xb = -0.1684;
    return 56680. / 243 + 32 * PI / (3 * std::sqrt(3)) + 128 * Xb / 9 - 16 * a1 + 32 * b1;
}

double BXsDecay::r6(double az, double bz, double a1, double b1) {
    double Xb = -0.1684;
    return 5710. / 729 - 16 * PI / (9 * std::sqrt(3)) - 64 * Xb / 27 - 10 * a1 / 3 + 44 * b1 / 3 + 12 * az + 20 * bz;
}

double BXsDecay::r8() {
    return 44. / 9 - 8 * PI2 / 27;
}

complex_t BXsDecay::G(double t) {
    if (t < 4) {
        return -2 * std::pow(std::atan(std::sqrt(t / (4 - t))), 2);
    }
    double L = std::log((std::sqrt(t) + sqrt(t - 4)) / 2);
    return complex_t {-PI2 / 2 + 2 * std::pow(L, 2), -2 * PI * L};
}

double BXsDecay::phi_11(double phi_22) {
    return phi_22 / 36;
}

double BXsDecay::phi_12(double phi_22) {
    return -phi_22 / 3;
}

double BXsDecay::phi_17(double phi_27) {
    return -phi_27 / 6;
}

double BXsDecay::phi_18(double phi_27) {
    return phi_27 / 18;
}

double BXsDecay::phi_22(double delta, double z) {
    auto i1 = [this, z] (double t) {
        return (1 - z * t) * std::pow(std::abs(G(t) / t + 0.5), 2);
    };

    auto i2 = [this, z] (double t) {
        return std::pow((1 - z * t) * std::abs(G(t) / t + 0.5), 2);
    };

    double c = (1 - delta) / z;
    double I1 = integrate(i1, 0, c, 1e-4);
    double I2 = integrate(i2, c, 1 / z, 1e-4);

    return 16 * z * (delta * I1 + I2) / 27;
}

double BXsDecay::phi_27(double delta, double z) {
    auto i1 = [this] (double t) {
        return std::real(G(t) + 0.5 * t);
    };

    auto i2 = [this, z] (double t) {
        return (1 - z * t) * std::real(G(t) + 0.5 * t);
    };

    double c = (1 - delta) / z;
    double I1 = integrate(i1, 0, c, 1e-4);
    double I2 = integrate(i2, c, 1 / z, 1e-4);

    return -8 * z * z * (delta * I1 + I2) / 9;
}

double BXsDecay::phi_28(double phi_27) {
    return -phi_27 / 3;
}

double BXsDecay::phi_47(double delta) {
    double phi47A=PI/54.*(3.*sqrt(3.)-PI)+delta*delta*delta/81.-25./108.*delta*delta+5./54.*delta+2./9.*(delta*delta+2.*delta+3.)*pow(atan(sqrt((1.-delta)/(3.+delta))),2.)-1./3.*(delta*delta+4.*delta+3.)*sqrt((1.-delta)/(3.+delta))*atan(sqrt((1.-delta)/(3.+delta)));
	double phi47B=(34.*delta*delta+59.*delta-18.)/486.*delta*delta*log(delta)/(1.-delta)+(433.*delta*delta*delta+429.*delta*delta-720.*delta)/2916.;
    
	return phi47A+phi47B;
}

double BXsDecay::phi_48(double phi_47) {
    return -phi_47 / 3;
}

double BXsDecay::phi_77(double delta) {
    double ld = std::log(delta);
    return -2 * std::pow(ld, 2) / 3 - 7 * ld / 3 - 31. / 9 + 10 * delta / 3 + std::pow(delta, 2) / 3 - 2 * std::pow(delta, 3) / 9 + delta * (delta - 4) * ld / 3;
}

double BXsDecay::phi_78(double delta) {
    return 8 * (Li2(1 - delta) - PI2 / 6 - delta * std::log(delta) + 9 * delta / 4 - std::pow(delta, 2) / 4 + std::pow(delta, 3) / 12) / 9;
}

double BXsDecay::phi_88(double delta, double mb_1S, double ms) {
    double u = 1 - delta;
    double a = delta * (delta + 2) + 4 * std::log(u);
    double b = 4 * Li2(u) - 2 * PI2 / 3 - delta * (delta + 2) * std::log(delta) + 8 * std::log(u) - 2 * std::pow(delta, 3) / 3 + 3 * std::pow(delta, 2) + 7 * delta;
    return (-2 * std::log(mb_1S / ms) * a + b) / 27;
}

double BXsDecay::K_i7_1(int i, double ri, double Lb, double phi_i7) {
    return ri - gamma_0[i - 1][6] * Lb / 2 + 2 * phi_i7;
}

double BXsDecay::K_77_1(double Lb, double phi_77) {
    return -182. / 9 + 8 * PI2 / 9 - gamma_0[6][6] * Lb + 4 * phi_77;
}

double BXsDecay::K_ij(int i, int j, double phi_ij) {
    return 2 * (1 + kron(i, j)) * phi_ij;
}

double BXsDecay::r22(double z) {
    double Lz = std::log(z);
    double Lz2 = Lz * Lz;
    double Lz3 = Lz2 * Lz;
    double Lz4 = Lz2 * Lz2;
    return  67454./6561.-124./729.*PI2
            -4./1215.*(11280.-1520.*PI2-171.*pow(PI,4.)-5760.*ZETA3+6840.*Lz-1440.*PI2*Lz-2520.*ZETA3*Lz+120.*Lz2+100.*Lz3-30.*Lz4)*z
	        -64./243.*PI2*(43.-12.*std::log(2.)-3.*Lz)*pow(z,1.5)
            -2./1215.*(11475.-380.*PI2+96.*pow(PI,4.)+7200.*ZETA3-1110.*Lz-1560.*PI2*Lz+1440.*ZETA3*Lz+990.*Lz2+260.*Lz3-60.*Lz4)*z*z
	        +2240./243.*PI2*pow(z,2.5)
            -2./2187.*(62471.-2424.*PI2-33264.*ZETA3-19494.*Lz-504.*PI2*Lz-5184.*Lz2+2160.*Lz3)*z*z*z
            -2464./6075.*PI2*pow(z,3.5)
            +(-15103841./546750.+7912./3645.*PI2+2368./81.*ZETA3+147038./6075.*Lz+352./243.*PI2*Lz+88./243.*Lz2-512./243.*Lz3)*z*z*z*z;
}

double BXsDecay::r22_large_z(double z) {
    double lz = std::log(z);
    return 27650. / 6561 + lz * (112. / 243 + 8 * lz / 9);
}

double BXsDecay::h22(double z, double delta) {
    return 0.01370+0.3357*delta-0.08668*delta*delta+(0.3575+1.825*delta-0.3743*delta*delta)*sqrt(z)+(-2.306-5.8*delta-6.226*delta*delta)*z+(3.449-0.548*delta+17.27*delta*delta)*pow(z,1.5);
}

double BXsDecay::h27(double z, double delta)
{
    return -0.1755-1.455*delta+1.119*delta*delta+(0.7260-7.23*delta+5.977*delta*delta)*sqrt(z)+(13.79+113.7*delta-100.4*delta*delta)*z+(-145.1-307.1*delta+388.5*delta*delta)*pow(z,1.5)+(475.2+313.*delta-775.8*delta*delta)*z*z+(-509.7-126.1*delta+646.2*delta*delta)*pow(z,2.5);
}

double BXsDecay::h28(double z, double delta)
{
    return 0.02605+0.1679*delta-0.197*delta*delta+(-0.03801+0.6017*delta-0.7558*delta*delta)*sqrt(z)+(2.755-10.03*delta+11.27*delta*delta)*z+(-27.05+68.47*delta-72.51*delta*delta)*pow(z,1.5)+(85.87-289.3*delta+297.7*delta*delta)*z*z+(-91.53+399.8*delta-399.9*delta*delta)*pow(z,2.5);
}

double BXsDecay::h88(double delta, double mb_1S, double ms) {
    return 4./27.*(((1.+0.5*delta)*delta*log(delta)-6.*log(1.-delta)-2.*Li2(1.-delta)+PI2/3.-16./3.*delta-5./3.*delta*delta+delta*delta*delta/9.)*log(mb_1S/ms)-2.*Li3(delta)+(5.-2.*log(delta))*(Li2(1.-delta)-PI2/6.)-PI2/12.*delta*(2.+delta)+(0.5*delta+0.25*delta*delta-log(1.-delta))*pow(log(delta),2.)+(151./18.-PI2/3.)*log(1.-delta)+(-53./12.-19./12.*delta+2./9.*delta*delta)*delta*log(delta)+787./72.*delta+227./72.*delta*delta-41./72.*delta*delta*delta);
}

double BXsDecay::F2nf(double z) {
    if(fpeq(z, 0.)) return 0.;
    
    return -std::log(1.-z)*std::log(1.-z)/(1.-z)/2.-13./36.*std::log(1.-z)/(1.-z)+(-PI2/18.+85./72.)/(1.-z)
	+(z*z-3.)/6./(z-1.)*Li2(1.-z)+(z*z-3.)/6./(z-1.)*std::log(1.-z)*std::log(z)-(1.+z)*std::log(1.-z)*std::log(1.-z)/4.-(6.*z*z-25.*z-1.)*std::log(1.-z)/36.+std::log(1.-z)/z/2.-(1.+z)*PI2/36.+(-49.+38.*z*z-55.*z)/72.;
}


double BXsDecay::h77(double delta) {
    auto f = [this] (double z) {
        return F2nf(z);
    };

    return 4 * integrate(f, 0, 1 - delta, 1e-3);
}
double BXsDecay::g(double m_c_3gev, double m_b_1S) {
    double rho = std::pow(m_c_3gev / m_b_1S, 2);
    return 1 + rho * (-8 + rho * (-12 * std::log(rho) + rho * (8 - rho)));
}

double BXsDecay::C(double g, double m_c_3gev, double mu_G2, double rho_D3, double rho_LS3) {
    double delta_as = QCDHelper::alpha_s(4.6, "running") - 0.22;
    double delta_b = QCDHelper::mass_b_msbar() - 4.18;
    double delta_c = m_c_3gev - 1;
    return g * (0.849 - 0.92 * delta_as + 0.0596 * delta_b - 0.2237 * delta_c - 0.0167 * mu_G2 - 0.203 * rho_D3 + 0.004 * rho_LS3);
}

double BXsDecay::x5(double P0, double K77_rem) {
    // return 17.35471;
    return P0 * K77_rem;
}

double BXsDecay::P0() {
    auto wilson = get_wilsons();
    complex_t C70 = wilson->getR(WGroup::B, WCoef::C7, QCDOrder::LO);
    complex_t C70p = wilson->getR(WGroup::BPrime, WCoef::CP7, QCDOrder::LO);
    return std::pow(std::abs(C70), 2) + std::pow(std::abs(C70p), 2);
}

double BXsDecay::P11() {
    auto C7 = get_wilsons()->getSR(WGroup::B, WCoef::C7);
    return 2 * std::real(C7[QCDOrder::LO] * std::conj(C7[QCDOrder::NLO]));
}

double BXsDecay::P12() {
    auto C7 = get_wilsons()->getSR(WGroup::B, WCoef::C7);
    return std::pow(std::abs(C7[QCDOrder::NLO]), 2) + 2 * std::real(C7[QCDOrder::LO] * std::conj(C7[QCDOrder::NNLO]));
}

template<std::size_t size>
double BXsDecay::gen_P00(const std::vector<scalar_t>& flat_K, const std::array<std::pair<int, int>, size>& indices) {
    auto wilson = get_wilsons();
    auto C = wilson->getAR(WGroup::B, QCDOrder::LO);
    auto Cp = wilson->getAR(WGroup::BPrime, QCDOrder::LO); 

    auto K = fill_K(flat_K, indices);
    double P {0};
    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 8; j++){
            P += std::real(C[static_cast<WCoef>(i)] * K[i][j] * std::conj(C[static_cast<WCoef>(j)]));
            P += std::real(Cp[static_cast<WCoef>(i + 12)] * K[i][j] * std::conj(Cp[static_cast<WCoef>(j + 12)]));
        }
    }
    
    return P;
}

template<std::size_t size>
double BXsDecay::gen_P01(const std::vector<scalar_t>& flat_K, const std::array<std::pair<int, int>, size>& indices) {
    auto wilson = get_wilsons();
    auto C0 = wilson->getAR(WGroup::B, QCDOrder::LO);
    auto C1 = wilson->getAR(WGroup::B, QCDOrder::NLO);

    auto K = fill_K(flat_K, indices);
    double P {0};
    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 8; j++){
            P += std::real(C0[static_cast<WCoef>(i)] * K[i][j] * std::conj(C1[static_cast<WCoef>(j)]));
        }
    }
    
    return 2 * P;
}

double BXsDecay::P22_rem(double x1,
                     double x2,
                     double x3,
                     double x4,
                     double x5,
                     complex_t r21,
                     double r22,
                     double dr21_dlogz)
{
    complex_t r21_0 {-1666. / 243, -80 * PI / 81};
    double r22_0 {67454. / 6561 - 124 * PI2 / 729}; 

    return x1 * (std::pow(std::abs(r21), 2) - std::pow(std::abs(r21_0), 2))
            + x2 * (r22 - r22_0)
            + x3 * std::real(r21 - r21_0)
            + x4 * dr21_dlogz
            + x5;
}

double BXsDecay::P(double alpha_mub,
                   double p0,
                   double p11,
                   double p12,
                   double p21,
                   double p22,
                   double p32)
{
    double k = alpha_mub / (4 * PI);
    return p0 + k * ((p11 + p21) + k * (p12 + p22 + p32));
}

double BXsDecay::Kc(double eta) {
    double kc = 0;
    for (size_t i = 0; i < a_i.size(); i++) {
        kc += d_i[i] * std::pow(eta, a_i[i]);
    }
    return kc;
}

double BXsDecay::Kt(double eta) {
    double f = std::pow(eta, 2. / 23);
    double f2 = f * f;
    auto wilson = get_wilsons();
    double C7 = std::real(wilson->getM(WGroup::B, WCoef::C7, QCDOrder::LO));
    double C8 = std::real(wilson->getM(WGroup::B, WCoef::C8, QCDOrder::LO));
    return (C7 + 23. / 36) * f2 - 8 * (C8 + 1. / 3) * (f2 - f) / 3;
}

double BXsDecay::r(double m_b_1S) {
    return QCDHelper::mass_b_msbar() / m_b_1S;
    // return QCDHelper::msbar_mass(5, winfo.matching_scale, "running") / m_b_1S;
}

double BXsDecay::N_eta_factor(double eta) {
    return std::pow(eta, 6. / 23) + std::pow(eta, -12. / 23);
}

double BXsDecay::N(double Kc,
                   double Kt,
                   double r,
                   double eta_factor,
                   double lambda_2,
                   double mc)
{
    return -(Kc + r * Kt) * eta_factor * lambda_2 / (18 * mc * mc);
}

double BXsDecay::k_SL(double alpha_mub) {
    return 2 * alpha_mub * std::log(winfo.matching_scale / winfo.hadronic_scale) / PI;
}

double BXsDecay::C2_em(double eta) {
    return -190 * std::pow(eta, -35. / 23) / 8073 - 359 * std::pow(eta, -17. / 23) / 3105 + 4276 * std::pow(eta, -12. / 23) / 121095
            + 350531 * std::pow(eta, -9. / 23) / 1009125 + 2 * std::pow(eta, -7. / 23) / 4347 - 5956 * std::pow(eta, 6. / 23) / 15525
            + 38380 * std::pow(eta, 14. / 23) / 169533 - 748 * std::pow(eta, 16. / 23) / 8625;
}

double BXsDecay::C8_em(double eta) {
    return -32 * std::pow(eta, -9. / 23) / 575 + 32 * std::pow(eta, -7. / 23) / 1449 + 640 * std::pow(eta, 14. / 23) / 1449 - 704 * std::pow(eta, 16. / 23) / 1725;
}

complex_t BXsDecay::C7_em(double eta, double C8_em, double C2_em) {
    auto wilson = get_wilsons();
    complex_t C7 = wilson->getM(WGroup::B, WCoef::C7, QCDOrder::LO);
    complex_t C8 = wilson->getM(WGroup::B, WCoef::C8, QCDOrder::LO);
    return (32 * std::pow(eta, -9. / 23) / 75 - 40 * std::pow(eta, -7. / 23) / 69 + 88 * std::pow(eta, 16. / 23) / 575) * C7 + C8_em * C8 + C2_em;
}

double BXsDecay::epsilon_em(double inv_alpha_em, double alpha_mub, double C7_em, double k) {
    complex_t C7 = get_wilsons()->getR(WGroup::B, WCoef::C7, QCDOrder::LO);
    return (2 * std::real(C7_em * std::conj(C7)) - k * std::pow(std::abs(C7), 2)) / (inv_alpha_em * alpha_mub);
}

double BXsDecay::ckm(double V_tb_r,
                     double V_tb_i,
                     double V_ts_r,
                     double V_ts_i,
                     double V_cb_r,
                     double V_cb_i)
{
    complex_t V_tb {V_tb_r, V_tb_i};
    complex_t V_ts_star {V_ts_r, -V_ts_i};
    complex_t V_cb {V_cb_r, V_cb_i};
    return std::pow(std::abs(V_ts_star * V_tb / V_cb), 2);
}

double BXsDecay::BR_B_Xs_gamma(double br_B__Xc_e_nu,
                               double ckm,
                               double inv_alpha_em,
                               double C,
                               double P,
                               double N,
                               double eps_em)
{
    return br_B__Xc_e_nu * ckm * 6 * (P + N + eps_em) / (inv_alpha_em * PI * C);
}

double BXsDecay::phi_11_b0(double phi_22_b0) {
    return phi_22_b0 / 36;
}

double BXsDecay::phi_12_b0(double phi_22_b0) {
    return -phi_22_b0 / 3;
}

double BXsDecay::phi_18_b0(double phi_28_b0) {
    return -phi_28_b0 / 6;
}

double BXsDecay::phi_ij_b0(double phi_ij_1, double b0, double Lb, double hij) {
    return b0 * (phi_ij_1 * Lb + hij);
}

double BXsDecay::K_17_b0(double K_27_b0) {
    return -K_27_b0 / 6;
}

double BXsDecay::K_27_b0(double b0, double r2, double az, double bz, double Lb, double phi_27_b0) {
    return b0 * std::real(-1.5 * r2 + 2 * (az + bz - 290. / 81) * Lb - 100 * Lb * Lb / 81) + 2 * phi_27_b0;
}

double BXsDecay::K_77_b0(double b0, double Lb, double phi_77_b0) {
    return b0 * (-3803. / 54 - 46 * PI2 / 27 + 80 * ZETA3 / 3 + (8 * PI2 / 9 - 98. / 3) * Lb - 16 * Lb * Lb / 3) + 4 * phi_77_b0;
}

double BXsDecay::K_78_b0(double b0, double Lb) {
    return b0 * (1256. / 81 - 64 * PI2 / 81 - 32 * ZETA3 / 9 + (-8 * PI2 / 27 + 188. / 27) * Lb + 8 * Lb * Lb / 9);
}

double BXsDecay::F2a(double z) {
    if(fpeq(z, 0.)) return 0;

    // double z2 = z * z;
    // double z3 = z2 * z;
    // double z4 = z3 * z;
    // double z5 = z4 * z;
    // double z6 = z5 * z;
    // double z7 = z6 * z;
    // double z8 = z7 * z;
    // double w = 1 - z;
    // double l1z = std::log(w);
    // double l1z2 = l1z * l1z;

	return 0.5*pow(log(1.-z),3.)/(1.-z)+21./8.*pow(log(1.-z),2.)/(1.-z)
	+(-PI2/6.+271./48.)*log(1.-z)/(1.-z)
	+(425./96.-PI2/6.-ZETA3/2.)/(1.-z)
	+(4.*z-4.*z*z+1.+z*z*z)/2./(z-1.)*(Li3(z/(2.-z))-Li3(-z/(2.-z))-2.*Li3(1./(2.-z))+ZETA3/4.)
	+((z*z*z-2.*z*z+2.*z-3.)/2./(z-1.)*log(1.-z)-(-140.*pow(z,4.)+219.*pow(z,3.)-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/12./z/pow(z-1.,3.))*Li2(z-1.)
-2.*pow(z-1.,2.)*Li3(z-1.)+((2.*pow(z,3.)-9.*z*z-2.*z+11.)/4./(z-1.)*log(1.-z)-(-27.*z*z+8.*pow(z,6.)-9.+21.*z-3.*z*z*z+64.*pow(z,4.)-46.*pow(z,5.))/12./z/pow(z-1.,3.))*Li2(1.-z)
-(-17.*z*z+4.*z+4.*z*z*z+11.)/4./(z-1.)*Li3(1.-z)-(2.*z*z*z+13.-9.*z*z)/4./(z-1.)*Li3(z)+(4.*z-4.*z*z+1.+z*z*z)/6./(z-1.)*pow(log(2.-z),3.) 
+(-(4.*z-4.*z*z+1.+z*z*z)/2./(z-1.)*pow(log(1.-z),2.)-(-140.*pow(z,4.)+219.*pow(z,3.)-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/12./z/pow(z-1.,3.)*log(1.-z) - (4.*z-4.*z*z+1.+z*z*z)/(z-1.)*PI2/12.)*log(2.-z)
+(z*z*z-2.*z*z+2.*z+1.)/4./z*pow(log(1.-z),3.)+(pow(z,5.)-3.*pow(z,4.)+5.*pow(z,3.)+7.*z*z+5.*z-9.)/24./z*pow(log(1.-z),2.)
+(-(z*z+8.*z-11.)/8./(z-1.)*pow(log(1.-z),2.)-(-27.*z*z+8.*pow(z,6.)-9.+21.*z-3.*z*z*z+64.*pow(z,4.)-46.*pow(z,5.))/12./z/pow(z-1.,3.)*log(1.-z))*log(z)
+((-z*z+z-3.)*PI2/12.-(4.*pow(z,5.)+151.*z+2.*pow(z,4.)-48.*z*z-41.*z*z*z-36.)/48./z*(z-1.))*log(1.-z)
-(z-2.)*(pow(z,4.)-z*z*z-11.*z*z+13.*z+3.)/z*PI2/72. + (z*z*z-11.*z*z-2.*z+18.)/4./(z-1.)*ZETA3-(8.*pow(z,4.)-244.*z*z*z+175.*z*z+598.*z-569.)/96./(z-1.);
}

double BXsDecay::F2na(double z) {
    if(fpeq(z, 0.)) return 0;

    // double z2 = z * z;
    // double z3 = z2 * z;
    // double z4 = z3 * z;
    // double z5 = z4 * z;
    // double z6 = z5 * z;
    // double z7 = z6 * z;
    // double z8 = z7 * z;
    // double w = 1 - z;
    // double l1z = std::log(w);
    // double l1z2 = l1z * l1z;

	return 11./8.*pow(log(1.-z),2.)/(1.-z)+(PI2/12.+95./144.)*log(1.-z)/(1.-z)+(ZETA3/4.-905./288.+17.*PI2/72.)/(1.-z)
	-(4.*z-4.*z*z+1.+z*z*z)/4./(z-1.)*(Li3(z/(2.-z))-Li3(-z/(2.-z))-2.*Li3(1./(2.-z))+ZETA3/4.)+pow(z-1.,2.)*Li3(z-1.)
	+(-(z*z*z-2.*z*z+2.*z-3.)/4./(z-1.)*log(1.-z)+(-140.*pow(z,4.)+219.*z*z*z-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/24./z/pow(z-1.,3.))*Li2(z-1.)
	+(z*(3.-z)/4.*log(1.-z)+(1.+z)*(2.*pow(z,4.)-29.*pow(z,3.)+73.*z*z-57.*z+15.)/24./pow(z-1.,3.))*Li2(1.-z)+(4.*z-4.*z*z+1.+z*z*z)/4./(z-1.)*Li3(z)
	+(z-3.)*z/2.*Li3(1.-z)-(4.*z-4.*z*z+1.+z*z*z)/12./(z-1.)*pow(log(2.-z),3.)+((4.*z-4.*z*z+1.+z*z*z)/4./(z-1.)*pow(log(1.-z),2.)
	+(-140.*pow(z,4.)+219.*pow(z,3.)-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/24./z/pow(z-1.,3.)*log(1.-z)+(4.*z-4.*z*z+1.+z*z*z)*PI2/24./(z-1.))*log(2.-z)
	+(1.+z)*(2.*pow(z,4.)-29.*z*z*z+73.*z*z-57.*z+15.)/24./pow(z-1.,3.)*log(1.-z)*log(z)-pow(z-1.,2.)/8.*pow(log(1.-z),3.)-(z+2.)*(z*z*z-5.*z*z+9.*z-35.)/48.*pow(log(1.-z),2.)
	+((z*z-z+3.)*PI2/24.+(6.*pow(z,5.)+72.-392.*pow(z,3.)+51.*pow(z,4.)+219.*z*z+92.*z)/144./z/(z-1.))*log(1.-z)
	+(pow(z,5.)-3.*pow(z,4.)-3.*pow(z,3.)+34.*z*z-24.*z+3.)/z*PI2/144.
	-(z*z*z-10.*z*z+6.*z+7.)/8./(z-1.)*ZETA3+(12.*pow(z,4.)-754.*pow(z,3.)+1191.*z*z+264.*z-761.)/288./(z-1.);
}

double BXsDecay::phi_77_rem_int(double delta) {
    auto f = [this] (double z) {
        return (16 * F2a(z) + 36 * F2na(z) + 87 * F2nf(z)) / 9.; 
    };

    return integrate(f, 0, 1 - delta, 1e-3);
}

double BXsDecay::phi_77_rem(double phi_77_int, double delta, double alpha_upsilon) {
    double ld = std::log(delta);
    return -4 * phi_77_int - 8 * PI * alpha_upsilon * (2 * delta * ld * ld + (4 + delta * (7 + delta * (-2 + delta))) * ld + 7 + delta * (-8. / 3 + delta * (-7 + delta * (4 - 4 * delta / 3)))) / (27 * delta);
}

double BXsDecay::K_11_rem(double K_22_rem) {
    return K_22_rem / 36;
}

double BXsDecay::K_12_rem(double K_22_rem) {
    return -K_22_rem / 6;
}

double BXsDecay::K_17_rem(double K_27_rem, double K_78_1, double LD) {
    return -K_27_rem / 6 + (5. / 16 - 3 * LD / 4) * K_78_1 - 1237. / 729 + 232 * ZETA3 / 27 + LD * (-20 + 70 * LD) / 27;
}

double BXsDecay::K_18_rem(double K_28_rem, double K_88_1, double LD) {
    return -K_28_rem / 6 + (5. / 16 - 3 * LD / 4) * K_88_1;
}

double BXsDecay::K_27_large_z(double LD) {
    return (218. / 243 - 208 * LD / 81);
}

double BXsDecay::K_22_rem(double K_27_1) {
    return std::pow(K_27_1, 2);
}

double BXsDecay::K_27_rem(double K_27_1, double K_77_1, double K_78_1, double K_47_1, double b0, double LD, double Lc, double Lb) {
    return K_27_1 * K_77_1 + (127. / 324 - 35 * LD / 27) * K_78_1 + 2 * (1 - LD) * (K_47_1 - b0 * (26. / 81 - 4 * Lb / 27)) / 3 + LD * (1150 - 4736 * LD) / 729 - 1617980. / 19683 + 20060 * ZETA3 / 243 + 1664 * Lc / 81;
}

double BXsDecay::K_28_rem(double K_27_1, double K_78_1, double K_88_1, double K_48_1, double LD) {
    return K_27_1 * K_78_1 + (127. / 324 - 35 * LD / 27) * K_88_1 + 2 * (1 - LD) * K_48_1 / 3;
}

double BXsDecay::K_77_rem(double K_77_1, double phi_77_1, double phi_77_rem, double z, double LD, double Lb, double alpha_upsilon) {
    return (K_77_1 - 4 * phi_77_1 + 2 * std::log(z) / 3) * K_77_1 + LD * (224. / 27 - 32 * LD / 9) - 79.2838955662 + Lb * (256 * PI2 / 27 - 2720. / 9 - 160 * Lb / 3) + 512 * PI * alpha_upsilon / 27 + 4 * phi_77_rem;
}

double BXsDecay::K_78_rem(double K_78_1, double LD) {
    return (-50. / 3 + 8 * PI2 / 3 - 2 * LD / 3) * K_78_1 + LD * (-112. / 81 + 16 * LD / 27) + 364. / 243;
}

double BXsDecay::K_88_rem(double K_88_1, double LD) {
    return (-50. / 3 + 8 * PI2 / 3 - 2 * LD / 3) * K_88_1;
}

double BXsDecay::c(double z) {
    return r2_large_z(z) + 1666. / 243;
}

double BXsDecay::target(double p_22, double x1, double x2, double x5, double z) {
    double a = std::pow(r2_large_z(z), 2) - std::pow(1666. / 243, 2) - std::pow(80 * PI / 81, 2);
    double b = std::real(r22_large_z(z) - 67454. / 6561 + 124 * PI2 / 729);
    return p_22 - x1 * a - x2 * b - x5;
}

double BXsDecay::x1() {
    auto wilson = get_wilsons();
    auto C = wilson->getAR(WGroup::B, QCDOrder::LO);
    auto CP = wilson->getAR(WGroup::BPrime, QCDOrder::LO);

    auto f = [] (complex_t c1, complex_t c2) { 
        return std::pow(std::abs(c1), 2) / 36 + std::pow(std::abs(c2), 2) - std::real(c1 * std::conj(c2)) / 3; 
    };

    return f(C[WCoef::C1], C[WCoef::C2]) + f(CP[WCoef::CP1], CP[WCoef::CP2]);
}

double BXsDecay::x2() {
    auto wilson = get_wilsons();
    auto C = wilson->getAllRunCoefficients(WGroup::B, QCDOrder::LO);
    auto CP = wilson->getAllRunCoefficients(WGroup::BPrime, QCDOrder::LO);

    auto f = [] (complex_t c1, complex_t c2, complex_t c7, complex_t c8) { 
        return std::real(c7 * std::conj(4019. * c1 / 486. - 1184. * c2 / 81. - 4. * c7 + 4. * c8 / 3.)); 
    };

    return f(C[WCoef::C1], C[WCoef::C2], C[WCoef::C7], C[WCoef::C8]) 
            + f(CP[WCoef::CP1], CP[WCoef::CP2], CP[WCoef::CP7], CP[WCoef::CP8]);
}

double BXsDecay::x3(double P22_z0, double P22_z1, double c_z0, double c_z1) {
    return (P22_z0 - P22_z1) / (c_z0 - c_z1);
}

double BXsDecay::x4(double P22_z0, double P22_z1, double c_z0, double c_z1) {
    double d = 208. / 81;
    return (c_z1 * P22_z0 - c_z0 * P22_z1) / (d * (c_z1 - c_z0));
}

void BXsDecay::build_op_tree() {

    // Formfactors and decay-specific parameters
    auto E0             = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 1));
    auto BR_B__Xc_e_nu  = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 2));
    auto mu_G2          = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 3));
    auto rho_D3         = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 4));
    auto rho_LS3        = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 5));
    auto lambda_2       = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 6));
    auto mu_c           = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 7));
    auto z0             = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 8));
    auto z1             = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Xs", 9));

    // SM parameters
    auto inv_alpha_em   = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 1));
    auto alpha_s_MZ     = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 3));
    auto M_Z            = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 4));
    auto mb_mb          = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 5));
    auto mt_pole        = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 6));
    auto V_cb_r         = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 12));
    auto V_cb_i         = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 12));
    auto V_tb_r         = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 22));
    auto V_tb_i         = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 22));
    auto V_ts_r         = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 21));
    auto V_ts_i         = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 21));
    auto m_d            = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 1));
    auto m_u            = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 2));
    auto m_s            = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 3));
    auto m_c            = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 4));

    auto dummy          = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));

    auto qcd            = std::make_shared<OperatorNode>("qcd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 0; });
    qcd->addChildren({alpha_s_MZ, M_Z, mt_pole, mb_mb, m_d, m_u, m_s, m_c});
    auto alpha_s_mu_b   = std::make_shared<OperatorNode>("alpha_s_mu_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->alpha_s_mub(); });
    alpha_s_mu_b->addChild(qcd);
    auto mb_1S          = std::make_shared<OperatorNode>("mb_1S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::mass_b_1S(); });
    mb_1S->addChild(qcd);
    auto eta            = std::make_shared<OperatorNode>("eta", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->eta(values[0]); });
    eta->addChild(alpha_s_mu_b);
    auto c2_em          = std::make_shared<OperatorNode>("C2_em", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->C2_em(values[0]); });
    c2_em->addChild(eta);
    auto c8_em          = std::make_shared<OperatorNode>("C8_em", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->C8_em(values[0]); });
    c8_em->addChild(eta);
    auto c7_em          = std::make_shared<OperatorNode>("C7_em", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->C7_em(values[0], values[1], values[2]); });
    c7_em->addChildren({eta, c8_em, c2_em, dummy});
    auto k_sl           = std::make_shared<OperatorNode>("k_SL", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->k_SL(values[0]); });
    k_sl->addChildren({alpha_s_mu_b});
    auto eps_em         = std::make_shared<OperatorNode>("eps_em", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->epsilon_em(values[0], values[1], values[2], values[3]); });
    eps_em->addChildren({inv_alpha_em, alpha_s_mu_b, c7_em, k_sl, dummy}); 
    auto k_c            = std::make_shared<OperatorNode>("K_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->Kc(values[0]); });
    k_c->addChildren({eta}); 
    auto k_t            = std::make_shared<OperatorNode>("K_t", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->Kt(values[0]); });
    k_t->addChildren({eta, dummy}); 
    auto f_N            = std::make_shared<OperatorNode>("f_N", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->N_eta_factor(values[0]); });
    f_N->addChildren({eta}); 
    auto r              = std::make_shared<OperatorNode>("r", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r(values[0]); });
    r->addChildren({mb_1S}); 
    auto N              = std::make_shared<OperatorNode>("N", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->N(values[0], values[1], values[2], values[3], values[4], values[5]); });
    N->addChildren({k_c, k_t, r, f_N, lambda_2, m_c}); 

    auto p_0            = std::make_shared<OperatorNode>("P_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->P0(); });
    p_0->addChildren({dummy});
    auto p_11           = std::make_shared<OperatorNode>("P_11", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->P11(); });
    p_11->addChildren({dummy});
    auto p_12           = std::make_shared<OperatorNode>("P_12", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->P12(); });
    p_12->addChildren({dummy});
    auto lb             = std::make_shared<OperatorNode>("L_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->Lb(values[0]); });
    lb->addChildren({mb_1S}); 
    auto mc_muc         = std::make_shared<OperatorNode>("m_c(mu_c)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->mc_muc(values[0]); });
    mc_muc->addChildren({mu_c}); 
    auto z              = std::make_shared<OperatorNode>("z", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->z(values[0], values[1]); });
    z->addChildren({mc_muc, mb_1S}); 
    auto delta          = std::make_shared<OperatorNode>("delta", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->delta(values[0], values[1]); });
    delta->addChildren({E0, mb_1S}); 
    auto az             = std::make_shared<OperatorNode>("a(z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->a(values[0]); });
    az->addChildren({z}); 
    auto bz             = std::make_shared<OperatorNode>("b(z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->b(values[0]); });
    bz->addChildren({z}); 
    auto r1             = std::make_shared<OperatorNode>("r_1(z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r1(values[0], values[1]); });
    r1->addChildren({az, bz});
    auto r2             = std::make_shared<OperatorNode>("r_2(z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r2(values[0], values[1]); });
    r2->addChildren({az, bz}); 
    auto r3             = std::make_shared<OperatorNode>("r_3", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r3(a(1), b(1)); });
    auto r4             = std::make_shared<OperatorNode>("r_4(z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r4(values[0], a(1), b(1)); });
    r4->addChildren({bz}); 
    auto r5             = std::make_shared<OperatorNode>("r_5", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r5(a(1), b(1)); }); 
    auto r6             = std::make_shared<OperatorNode>("r_6(z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r6(values[0], values[1], a(1), b(1)); });
    r6->addChildren({az, bz}); 
    auto r8             = std::make_shared<OperatorNode>("r_8", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r8(); }); 
    auto phi_22         = std::make_shared<OperatorNode>("phi_22", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_22(values[0], values[1]); });
    phi_22->addChildren({delta, z}); 
    auto phi_27         = std::make_shared<OperatorNode>("phi_27", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_27(values[0], values[1]); });
    phi_27->addChildren({delta, z}); 
    auto phi_47         = std::make_shared<OperatorNode>("phi_47", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_47(values[0]); });
    phi_47->addChildren({delta}); 
    auto phi_77         = std::make_shared<OperatorNode>("phi_77", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_77(values[0]); });
    phi_77->addChildren({delta});
    auto phi_78         = std::make_shared<OperatorNode>("phi_78", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_78(values[0]); });
    phi_78->addChildren({delta});
    auto phi_88         = std::make_shared<OperatorNode>("phi_88", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_88(values[0], values[1], values[2]); });
    phi_88->addChildren({delta, mb_1S, m_s});
    auto phi_11         = std::make_shared<OperatorNode>("phi_11", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_11(values[0]); });
    phi_11->addChildren({phi_22});
    auto phi_12         = std::make_shared<OperatorNode>("phi_12", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_12(values[0]); });
    phi_12->addChildren({phi_22});
    auto phi_17         = std::make_shared<OperatorNode>("phi_17", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_17(values[0]); });
    phi_17->addChildren({phi_27});
    auto phi_18         = std::make_shared<OperatorNode>("phi_18", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_18(values[0]); });
    phi_18->addChildren({phi_27});
    auto phi_28         = std::make_shared<OperatorNode>("phi_28", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_28(values[0]); });
    phi_28->addChildren({phi_27});
    auto phi_48         = std::make_shared<OperatorNode>("phi_48", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_48(values[0]); });
    phi_48->addChildren({phi_47});
    auto k_11           = std::make_shared<OperatorNode>("K_11", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(1, 1, values[0]); });
    k_11->addChildren({phi_11});
    auto k_12           = std::make_shared<OperatorNode>("K_12", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(1, 2, values[0]); });
    k_12->addChildren({phi_12});
    auto k_17           = std::make_shared<OperatorNode>("K_17", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(1, values[0], values[1], values[2]); });
    k_17->addChildren({r1, lb, phi_17});
    auto k_18           = std::make_shared<OperatorNode>("K_18", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(1, 8, values[0]); });
    k_18->addChildren({phi_18});
    auto k_22           = std::make_shared<OperatorNode>("K_22", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(2, 2, values[0]); });
    k_22->addChildren({phi_22});
    auto k_27           = std::make_shared<OperatorNode>("K_27", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(2, values[0], values[1], values[2]); });
    k_27->addChildren({r2, lb, phi_27});
    auto k_28           = std::make_shared<OperatorNode>("K_28", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(2, 8, values[0]); });
    k_28->addChildren({phi_28});
    auto k_37           = std::make_shared<OperatorNode>("K_37", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(3, values[0], values[1], 0); });
    k_37->addChildren({r3, lb});
    auto k_47           = std::make_shared<OperatorNode>("K_47", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(4, values[0], values[1], values[2]); });
    k_47->addChildren({r4, lb, phi_47});
    auto k_48           = std::make_shared<OperatorNode>("K_48", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(4, 8, values[0]); });
    k_48->addChildren({phi_48});
    auto k_57           = std::make_shared<OperatorNode>("K_57", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(5, values[0], values[1], 0); });
    k_57->addChildren({r5, lb});
    auto k_67           = std::make_shared<OperatorNode>("K_67", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(6, values[0], values[1], 0); });
    k_67->addChildren({r6, lb});
    auto k_77           = std::make_shared<OperatorNode>("K_77", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_77_1(values[0], values[1]); });
    k_77->addChildren({lb, phi_77});
    auto k_78           = std::make_shared<OperatorNode>("K_78", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_i7_1(8, values[0], values[1], values[2]); });
    k_78->addChildren({r8, lb, phi_78});
    auto k_88           = std::make_shared<OperatorNode>("K_88", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(8, 8, values[0]); });
    k_88->addChildren({phi_88});
    auto p_21           = std::make_shared<OperatorNode>("P_21", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->gen_P00(values, nonzero_K1); });
    p_21->addChildren({k_11, k_12, k_17, k_18, k_22, k_27, k_28, k_37, k_47, k_48, k_57, k_67, k_77, k_78, k_88, dummy});
    auto p_32           = std::make_shared<OperatorNode>("P_32", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->gen_P01(values, nonzero_K1); });
    p_32->addChildren({k_11, k_12, k_17, k_18, k_22, k_27, k_28, k_37, k_47, k_48, k_57, k_67, k_77, k_78, k_88, dummy});
    auto h_22           = std::make_shared<OperatorNode>("h_22", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->h22(values[0], values[1]); });
    h_22->addChildren({z, delta});
    auto h_27           = std::make_shared<OperatorNode>("h_27", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->h27(values[0], values[1]); });
    h_27->addChildren({z, delta});
    auto h_28           = std::make_shared<OperatorNode>("h_28", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->h28(values[0], values[1]); });
    h_28->addChildren({z, delta});
    auto h_77           = std::make_shared<OperatorNode>("h_77", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->h77(values[0]); });
    h_77->addChildren({delta});
    auto h_88           = std::make_shared<OperatorNode>("h_88", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->h88(values[0], values[1], values[2]); });
    h_88->addChildren({delta, mb_1S, m_s});
    auto beta_0         = std::make_shared<OperatorNode>("beta_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::beta[QCDHelper::get_nf(winfo.hadronic_scale)][0]; });
    beta_0->addChild(qcd);
    auto phi_22_b0      = std::make_shared<OperatorNode>("phi_22_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_ij_b0(values[0], values[1], values[2], values[3]); });
    phi_22_b0->addChildren({phi_22, beta_0, lb, h_22});
    auto phi_11_b0      = std::make_shared<OperatorNode>("phi_11_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_11_b0(values[0]); });
    phi_11_b0->addChildren({phi_22_b0});
    auto phi_12_b0      = std::make_shared<OperatorNode>("phi_12_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_12_b0(values[0]); });
    phi_12_b0->addChildren({phi_22_b0});
    auto phi_27_b0      = std::make_shared<OperatorNode>("phi_27_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_ij_b0(values[0], values[1], values[2], values[3]); });
    phi_27_b0->addChildren({phi_27, beta_0, lb, h_27});
    auto phi_28_b0      = std::make_shared<OperatorNode>("phi_28_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_ij_b0(values[0], values[1], values[2], values[3]); });
    phi_28_b0->addChildren({phi_28, beta_0, lb, h_28});
    auto phi_18_b0      = std::make_shared<OperatorNode>("phi_18_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_18_b0(values[0]); });
    phi_18_b0->addChildren({phi_28_b0});
    auto phi_77_b0      = std::make_shared<OperatorNode>("phi_77_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_ij_b0(values[0], values[1], values[2], values[3]); });
    phi_77_b0->addChildren({phi_77, beta_0, lb, h_77});
    auto phi_88_b0      = std::make_shared<OperatorNode>("phi_88_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_ij_b0(values[0], values[1], values[2], values[3]); });
    phi_88_b0->addChildren({phi_88, beta_0, lb, h_88});
    auto r22            = std::make_shared<OperatorNode>("r_22", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->r22(values[0]); });
    r22->addChildren({z});
    auto k_11_b0        = std::make_shared<OperatorNode>("K_11_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(1, 1, values[0]); });
    k_11_b0->addChildren({phi_11_b0});
    auto k_12_b0        = std::make_shared<OperatorNode>("K_12_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(1, 2, values[0]); });
    k_12_b0->addChildren({phi_12_b0});
    auto k_27_b0        = std::make_shared<OperatorNode>("K_27_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_27_b0(values[0], values[1], values[2], values[3], values[4], values[5]); });
    k_27_b0->addChildren({beta_0, r22, az, bz, lb, phi_27_b0});
    auto k_17_b0        = std::make_shared<OperatorNode>("K_17_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_17_b0(values[0]); });
    k_17_b0->addChildren({k_27_b0});
    auto k_18_b0        = std::make_shared<OperatorNode>("K_18_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(1, 8, values[0]); });
    k_18_b0->addChildren({phi_18_b0});
    auto k_22_b0        = std::make_shared<OperatorNode>("K_22_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(2, 2, values[0]); });
    k_22_b0->addChildren({phi_22_b0});
    auto k_28_b0        = std::make_shared<OperatorNode>("K_28_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(2, 8, values[0]); });
    k_28_b0->addChildren({phi_28_b0});
    auto k_77_b0        = std::make_shared<OperatorNode>("K_77_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_77_b0(values[0], values[1], values[2]); });
    k_77_b0->addChildren({beta_0, lb, phi_77_b0});
    auto k_78_b0        = std::make_shared<OperatorNode>("K_78_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_78_b0(values[0], values[1]); });
    k_78_b0->addChildren({beta_0, lb});
    auto k_88_b0        = std::make_shared<OperatorNode>("K_88_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_ij(8, 8, values[0]); });
    k_88_b0->addChildren({phi_88_b0});
    auto p_22_b0        = std::make_shared<OperatorNode>("P_22_b0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->gen_P00(values, nonzero_K2); });
    p_22_b0->addChildren({k_11_b0, k_12_b0, k_17_b0, k_18_b0, k_22_b0, k_27_b0, k_28_b0, k_77_b0, k_78_b0, k_88_b0, dummy}); 

    auto x_1            = std::make_shared<OperatorNode>("x_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->x1(); });
    x_1->addChildren({dummy});
    auto x_2            = std::make_shared<OperatorNode>("x_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->x2(); });
    x_2->addChildren({dummy});
    auto ld             = std::make_shared<OperatorNode>("LD", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->LD(values[0], values[1]); });
    ld->addChildren({lb, z});
    auto ld_z0          = std::make_shared<OperatorNode>("LD_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->LD(values[0], values[1]); });
    ld_z0->addChildren({lb, z0});
    auto ld_z1          = std::make_shared<OperatorNode>("LD_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->LD(values[0], values[1]); });
    ld_z1->addChildren({lb, z1});
    auto lc             = std::make_shared<OperatorNode>("Lc", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->Lc(values[0], values[1]); });
    lc->addChildren({mu_c, mc_muc});
    auto alpha_s_ups    = std::make_shared<OperatorNode>("alpha_s_upsilon", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->alpha_s_upsilon(values[0]); });
    alpha_s_ups->addChild(mb_1S);
    auto I_77           = std::make_shared<OperatorNode>("I_77", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_77_rem_int(values[0]); });
    I_77->addChild(delta);
    auto phi_77_rem     = std::make_shared<OperatorNode>("phi_77_rem", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->phi_77_rem(values[0], values[1], values[2]); });
    phi_77_rem->addChildren({I_77, delta, alpha_s_ups});
    auto k77_rem_z      = std::make_shared<OperatorNode>("k_77_rem_z", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_77_rem(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    k77_rem_z->addChildren({k_77, phi_77, phi_77_rem, z, ld, lb, alpha_s_ups});
    auto x_5            = std::make_shared<OperatorNode>("x_5", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->x5(values[0], values[1]); });
    x_5->addChildren({p_0, k77_rem_z});

    auto k_27_z0         = std::make_shared<OperatorNode>("K_27_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_27_large_z(values[0]); });
    k_27_z0->addChildren({ld_z0});        
    auto k_22_rem_z0     = std::make_shared<OperatorNode>("K_22_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_22_rem(values[0]); });
    k_22_rem_z0->addChildren({k_27_z0});
    auto k_11_rem_z0     = std::make_shared<OperatorNode>("K_11_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_11_rem(values[0]); });
    k_11_rem_z0->addChildren({k_22_rem_z0});
    auto k_12_rem_z0     = std::make_shared<OperatorNode>("K_12_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_12_rem(values[0]); });
    k_12_rem_z0->addChildren({k_22_rem_z0});
    auto k_27_rem_z0     = std::make_shared<OperatorNode>("K_27_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_27_rem(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]); });
    k_27_rem_z0->addChildren({k_27_z0, k_77, k_78, k_47, beta_0, ld_z0, lc, lb});
    auto k_28_rem_z0     = std::make_shared<OperatorNode>("K_28_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_28_rem(values[0], values[1], values[2], values[3], values[4]); });
    k_28_rem_z0->addChildren({k_27_z0, k_78, k_88, k_48, ld_z0});
    auto k_18_rem_z0     = std::make_shared<OperatorNode>("K_18_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_18_rem(values[0], values[1], values[2]); });
    k_18_rem_z0->addChildren({k_28_rem_z0, k_88, ld_z0});
    auto k_17_rem_z0     = std::make_shared<OperatorNode>("K_17_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_17_rem(values[0], values[1], values[2]); });
    k_17_rem_z0->addChildren({k_27_rem_z0, k_78, ld_z0});
    auto k_77_rem_z0     = std::make_shared<OperatorNode>("K_77_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_77_rem(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    k_77_rem_z0->addChildren({k_77, phi_77, phi_77_rem, z0, ld_z0, lb, alpha_s_ups});
    auto k_78_rem_z0     = std::make_shared<OperatorNode>("K_78_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_78_rem(values[0], values[1]); });
    k_78_rem_z0->addChildren({k_78, ld_z0});
    auto k_88_rem_z0     = std::make_shared<OperatorNode>("K_88_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_88_rem(values[0], values[1]); });
    k_88_rem_z0->addChildren({k_88, ld_z0});
    auto p_22_rem_z0     = std::make_shared<OperatorNode>("p_22_rem_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->gen_P00(values, nonzero_K2); });
    p_22_rem_z0->addChildren({k_11_rem_z0, k_12_rem_z0, k_17_rem_z0, k_18_rem_z0, k_22_rem_z0, k_27_rem_z0, k_28_rem_z0, k_77_rem_z0, k_78_rem_z0, k_88_rem_z0, dummy});
    auto target_z0      = std::make_shared<OperatorNode>("target_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->target(values[0], values[1], values[2], values[3], values[4]); });
    target_z0->addChildren({p_22_rem_z0, x_1, x_2, x_5, z0});

    auto k_27_z1         = std::make_shared<OperatorNode>("K_27_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_27_large_z(values[0]); });
    k_27_z1->addChildren({ld_z1});   
    auto k_22_rem_z1     = std::make_shared<OperatorNode>("K_22_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_22_rem(values[0]); });
    k_22_rem_z1->addChildren({k_27_z1});
    auto k_11_rem_z1     = std::make_shared<OperatorNode>("K_11_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_11_rem(values[0]); });
    k_11_rem_z1->addChildren({k_22_rem_z1});
    auto k_12_rem_z1     = std::make_shared<OperatorNode>("K_12_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_12_rem(values[0]); });
    k_12_rem_z1->addChildren({k_22_rem_z1});
    auto k_27_rem_z1     = std::make_shared<OperatorNode>("K_27_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_27_rem(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]); });
    k_27_rem_z1->addChildren({k_27_z1, k_77, k_78, k_47, beta_0, ld_z1, lc, lb});
    auto k_28_rem_z1     = std::make_shared<OperatorNode>("K_28_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_28_rem(values[0], values[1], values[2], values[3], values[4]); });
    k_28_rem_z1->addChildren({k_27_z1, k_78, k_88, k_48, ld_z1});
    auto k_18_rem_z1     = std::make_shared<OperatorNode>("K_18_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_18_rem(values[0], values[1], values[2]); });
    k_18_rem_z1->addChildren({k_28_rem_z1, k_88, ld_z1});
    auto k_17_rem_z1     = std::make_shared<OperatorNode>("K_17_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_17_rem(values[0], values[1], values[2]); });
    k_17_rem_z1->addChildren({k_27_rem_z1, k_78, ld_z1});
    auto k_77_rem_z1     = std::make_shared<OperatorNode>("K_77_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_77_rem(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    k_77_rem_z1->addChildren({k_77, phi_77, phi_77_rem, z1, ld_z1, lb, alpha_s_ups});
    auto k_78_rem_z1     = std::make_shared<OperatorNode>("K_77_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_78_rem(values[0], values[1]); });
    k_78_rem_z1->addChildren({k_78, ld_z1});
    auto k_88_rem_z1     = std::make_shared<OperatorNode>("K_88_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->K_88_rem(values[0], values[1]); });
    k_88_rem_z1->addChildren({k_88, ld_z1});
    auto p_22_rem_z1     = std::make_shared<OperatorNode>("p_22_rem_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->gen_P00(values, nonzero_K2); });
    p_22_rem_z1->addChildren({k_11_rem_z1, k_12_rem_z1, k_17_rem_z1, k_18_rem_z1, k_22_rem_z1, k_27_rem_z1, k_28_rem_z1, k_77_rem_z1, k_78_rem_z1, k_88_rem_z1, dummy});
    auto target_z1      = std::make_shared<OperatorNode>("target_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->target(values[0], values[1], values[2], values[3], values[4]); });
    target_z1->addChildren({p_22_rem_z1, x_1, x_2, x_5, z1});

    auto c_z0           = std::make_shared<OperatorNode>("c_z0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->c(values[0]); });
    c_z0->addChild(z0);
    auto c_z1           = std::make_shared<OperatorNode>("c_z1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->c(values[0]); });
    c_z1->addChild(z1);
    auto x_3            = std::make_shared<OperatorNode>("x_3", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->x3(values[0], values[1], values[2], values[3]); });
    x_3->addChildren({target_z0, target_z1, c_z0, c_z1});
    auto x_4            = std::make_shared<OperatorNode>("x_4", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->x4(values[0], values[1], values[2], values[3]); });
    x_4->addChildren({target_z0, target_z1, c_z0, c_z1});
    auto dr21_dlogz     = std::make_shared<OperatorNode>("d(r_21)/d(log z)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->dr2_dlogz(values[0]); });
    dr21_dlogz->addChildren({z});
    auto p_22_rem       = std::make_shared<OperatorNode>("P_22_rem", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->P22_rem(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]); });
    p_22_rem->addChildren({x_1, x_2, x_3, x_4, x_5, r2, r22, dr21_dlogz});
    auto p_22           = std::make_shared<OperatorNode>("P_22", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    p_22->addChildren({p_22_b0, p_22_rem});
    auto P              = std::make_shared<OperatorNode>("P", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->P(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    P->addChildren({alpha_s_mu_b, p_0, p_11, p_12, p_21, p_22, p_32});
    auto mc_3gev        = std::make_shared<OperatorNode>("m_c(3 GeV)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->mc_3GeV(); });
    mc_3gev->addChildren({qcd}); 
    auto g              = std::make_shared<OperatorNode>("g(rho)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->g(values[0], values[1]); });
    g->addChildren({mc_3gev, mb_1S}); 
    auto C              = std::make_shared<OperatorNode>("C", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->C(values[0], values[1], values[2], values[3], values[4]); });
    C->addChildren({g, mc_3gev, mu_G2, rho_D3, rho_LS3, qcd}); 
    auto ckm            = std::make_shared<OperatorNode>("ckm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0], values[1], values[2], values[3], values[4], values[5]); });
    ckm->addChildren({V_tb_r, V_tb_i, V_ts_r, V_ts_i, V_cb_r, V_cb_i});
    auto br_B_Xs_gamma  = std::make_shared<OperatorNode>("BR(B > Xs gamma)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_B_Xs_gamma(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    br_B_Xs_gamma->addChildren({BR_B__Xc_e_nu, ckm, inv_alpha_em, C, P, N, eps_em});

    roots.emplace(Observables::BR_B_XS_GAMMA, br_B_Xs_gamma);
}