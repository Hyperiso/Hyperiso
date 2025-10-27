#include "BXsDecay.h"

void BXsDecay::load_params() {
    ObsParameterProxy p;
    cache.alpha_em = 1. / p(ParamId{ParameterType::SM, "SMINPUTS", 1});
    cache.m_s = p(ParamId{ParameterType::SM, "MASS", 3});
    cache.m_c = p(ParamId{ParameterType::SM, "MASS", 4});
    cache.m_b_mb = p(ParamId{ParameterType::SM, "QCD", {5, 1}});
    cache.m_b_1S = p(ParamId{ParameterType::SM, "QCD", {5, 3}});
    cache.ckm_factor = std::pow(std::abs(std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}) / p(ParamId{ParameterType::SM, "VCKM", {1, 2}})), 2);
    cache.mu_b = w_config.hadronic_scale;
    cache.mu_W = w_config.matching_scale;
    cache.beta_0 = ObsQCDProxy().get_constants()->beta[5][0]; // TODO : compute n_f based on scale ?
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.alpha_s_upsilon = ObsQCDProxy()(AlphasConfig(cache.m_b_1S, MassType::POLE, MassType::POLE));
    cache.eta = ObsQCDProxy()(AlphasConfig(cache.mu_W, MassType::POLE, MassType::POLE)) / cache.alpha_s_mu_b;
    cache.E0 = p(ParamId{ParameterType::DECAY, "B_Xs", 1});
    cache.BR_B__Xc_e_nu_exp = p(ParamId{ParameterType::DECAY, "B_Xs", 2});
    cache.mu_G2 = p(ParamId{ParameterType::DECAY, "B_Xs", 3});
    cache.rho_D3= p(ParamId{ParameterType::DECAY, "B_Xs", 4});
    cache.rho_LS3= p(ParamId{ParameterType::DECAY, "B_Xs", 5});
    cache.lambda_2= p(ParamId{ParameterType::DECAY, "B_Xs", 6});
    cache.mu_c= p(ParamId{ParameterType::DECAY, "B_Xs", 7});
    cache.z0= p(ParamId{ParameterType::DECAY, "B_Xs", 8});
    cache.z1= p(ParamId{ParameterType::DECAY, "B_Xs", 9});

    cache.r_msbar_1S = ObsQCDProxy()(MassConfig(5, cache.mu_W, MassType::MSBAR, MassType::POLE)) / cache.m_b_1S; // mu_W or mu_b ????
    cache.m_c_mu_c = ObsQCDProxy()(MassConfig(4, cache.mu_c, MassType::MSBAR, MassType::POLE));
    cache.m_c_3gev = ObsQCDProxy()(MassConfig(4, 3.0, MassType::MSBAR, MassType::POLE));
    cache.z = std::pow(cache.m_c_mu_c / cache.m_b_1S, 2);
    cache.delta = 1 - 2 * cache.E0 / cache.m_b_1S;
    cache.L_b = 2 * std::log(cache.mu_b / cache.m_b_1S);

    cache.L_c = 2 * std::log(cache.mu_c / cache.m_c_mu_c);
    
    cache.C_b_LO = w_proxy->getAR(WGroup::B, QCDOrder::LO);
    auto CP_b_LO = w_proxy->getAR(WGroup::BPrime, QCDOrder::LO);
    cache.C_b_LO.insert(CP_b_LO.begin(), CP_b_LO.end());
    cache.C_b_NLO = w_proxy->getAR(WGroup::B, QCDOrder::NLO);
    cache.C_b_NNLO = w_proxy->getAR(WGroup::B, QCDOrder::NNLO);
    cache.C_w = w_proxy->getAM(WGroup::B, QCDOrder::LO);
}   

double BXsDecay::gen_P00(const std::array<std::array<double, 8>, 8>& K) {
    double P {0};
    auto C_ids = WCoefMapper::get_group(WGroup::B);
    auto CP_ids = WCoefMapper::get_group(WGroup::BPrime);
    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 8; j++) {
            P += std::real(cache.C_b_LO[C_ids[i]] * K[i][j] * std::conj(cache.C_b_LO[C_ids[j]]));
            P += std::real(cache.C_b_LO[CP_ids[i]] * K[i][j] * std::conj(cache.C_b_LO[CP_ids[j]]));
        }
    }
    
    return P;
}

double BXsDecay::gen_P01(const std::array<std::array<double, 8>, 8>& K) {
    double P {0};
    auto C_ids = WCoefMapper::get_group(WGroup::B);
    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 8; j++) {
            P += std::real(cache.C_b_LO[C_ids[i]] * K[i][j] * std::conj(cache.C_b_NLO[C_ids[j]]));
        }
    }
    
    return 2 * P;
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

scalar_t BXsDecay::G(double t) {
    if (t < 4) {
        return -2 * std::pow(std::atan(std::sqrt(t / (4 - t))), 2);
    }
    double L = std::log((std::sqrt(t) + sqrt(t - 4)) / 2);
    return scalar_t {-PI2 / 2 + 2 * std::pow(L, 2), -2 * PI * L};
}

double BXsDecay::phi_22(double z, double delta) {
    auto i1 = [this, z] (double t) {
        return (1 - z * t) * std::pow(std::abs(G(t) / t + 0.5), 2);
    };

    auto i2 = [this, z] (double t) {
        return std::pow((1 - z * t) * std::abs(G(t) / t + 0.5), 2);
    };

    double c = (1 - delta) / z;
    double I1 = integrate(i1, 0, c, 1e-4);
    double I2 = integrate(i2, c, 1 / z, 1e-4);

    return 16. * z * (delta * I1 + I2) / 27.;
}

double BXsDecay::phi_27(double z, double delta) {
    auto i1 = [this] (double t) {
        return std::real(G(t) + 0.5 * t);
    };

    auto i2 = [this, z] (double t) {
        return (1 - z * t) * std::real(G(t) + 0.5 * t);
    };

    double c = (1 - delta) / z;
    double I1 = integrate(i1, 0, c, 1e-4);
    double I2 = integrate(i2, c, 1 / z, 1e-4);

    return -8. * std::pow(z, 2) * (delta * I1 + I2) / 9.;
}

double BXsDecay::phi_47(double delta) {
    double d = delta;
    double phi47A=PI/54.*(3.*sqrt(3.)-PI)+d*d*d/81.-25./108.*d*d+5./54.*d+2./9.*(d*d+2.*d+3.)*pow(atan(sqrt((1.-d)/(3.+d))),2.)-1./3.*(d*d+4.*d+3.)*sqrt((1.-d)/(3.+d))*atan(sqrt((1.-d)/(3.+d)));
	double phi47B=(34.*d*d+59.*d-18.)/486.*d*d*log(d)/(1.-d)+(433.*d*d*d+429.*d*d-720.*d)/2916.;
	return phi47A+phi47B;
}
double BXsDecay::phi_77(double delta) {
    double ld = std::log(delta);
    return -2 * std::pow(ld, 2) / 3 - 7 * ld / 3 - 31. / 9 + 10 * delta / 3 + std::pow(delta, 2) / 3 - 2 * std::pow(delta, 3) / 9 + delta * (delta - 4) * ld / 3;
}

double BXsDecay::phi_78(double delta) {
    return 8 * (Li2(1 - delta) - PI2 / 6 - delta * std::log(delta) + 9 * delta / 4 - std::pow(delta, 2) / 4 + std::pow(delta, 3) / 12) / 9;
}

double BXsDecay::phi_88(double delta) {
    double u = 1 - delta;
    double a = delta * (delta + 2) + 4 * std::log(u);
    double b = 4 * Li2(u) - 2 * PI2 / 3 - delta * (delta + 2) * std::log(delta) + 8 * std::log(u) - 2 * std::pow(delta, 3) / 3 + 3 * std::pow(delta, 2) + 7 * delta;
    return (-2 * std::log(cache.m_b_1S / cache.m_s) * a + b) / 27;
}

std::array<std::array<double, 8>, 8> BXsDecay::phi_1(double delta, double z) {
    std::array<std::array<double, 8>, 8> phi {};

    phi[1][1] = phi_22(z, delta);
    phi[1][6] = phi_27(z, delta);
    phi[3][6] = phi_47(delta);
    phi[6][6] = phi_77(delta);
    phi[6][7] = phi_78(delta);
    phi[7][7] = phi_88(delta);
    phi[0][0] = phi[1][1] / 36.0;
    phi[0][1] = -phi[1][1] / 3.0;
    phi[0][6] = -phi[1][6] / 6.0;
    phi[0][7] = phi[1][6] / 18.0;
    phi[1][7] = -phi[1][6] / 3.0;
    phi[3][7] = -phi[3][6] / 3.0;
    return phi;
}

std::array<double, 8> BXsDecay::r_1(double z) {
    std::array<double, 8> r {};
    double az = a(z);
    double bz = b(z);
    double a1 = a(1.0);
    double b1 = b(1.0);
    r[0] = 833. / 729 - (az + bz) / 3;
    r[1] = -1666. / 243 + 2 * (az + bz);
    r[2] = 2392. / 243 + 8 * PI / (3 * std::sqrt(3)) + 32 * cache.X_b / 9 - a1 + 2 * b1;
    r[3] = -761. / 729 - 4 * PI / (9 * std::sqrt(3)) - 16 * cache.X_b / 27 + a1 / 6 + 5 * b1 /3 + 2 * bz;
    r[4] = 56680. / 243 + 32 * PI / (3 * std::sqrt(3)) + 128 * cache.X_b / 9 - 16 * a1 + 32 * b1;
    r[5] = 5710. / 729 - 16 * PI / (9 * std::sqrt(3)) - 64 * cache.X_b / 27 - 10 * a1 / 3 + 44 * b1 / 3 + 12 * az + 20 * bz;
    r[6] = -182. / 9 + 8 * PI2 / 9;
    r[7] = 44. / 9 - 8 * PI2 / 27;
    return r;
}

std::array<std::array<double, 8>, 8> BXsDecay::K_1() {
    std::array<std::array<double, 8>, 8> K {};
    auto r = r_1(cache.z);
    auto phi = phi_1(cache.delta, cache.z);

    for (size_t i = 0; i < 8; i++) {
        for (size_t j = i; j < 8; j++) {
            K[i][j] = 2. * (1 + kron(i, j)) * phi[i][j];
        }
    }

    for (size_t i = 0; i < 6; i++)
        K[i][6] += r[i] - 0.5 * gamma_i7[i] * cache.L_b;

    K[6][6] += r[6] - gamma_i7[6] * cache.L_b;
    K[6][7] += r[7] - 0.5 * gamma_i7[7] * cache.L_b;

    for (size_t i = 0; i < 8; i++) {
        for (size_t j = i; j < 8; j++)  {
            K[j][i] = K[i][j];
        }
    }

    return K;
}

double BXsDecay::F2nf(double z) {
    if(fpeq(z, 0.)) return 0.;
    
    return -std::log(1.-z)*std::log(1.-z)/(1.-z)/2.-13./36.*std::log(1.-z)/(1.-z)+(-PI2/18.+85./72.)/(1.-z)
	+(z*z-3.)/6./(z-1.)*Li2(1.-z)+(z*z-3.)/6./(z-1.)*std::log(1.-z)*std::log(z)-(1.+z)*std::log(1.-z)*std::log(1.-z)/4.-(6.*z*z-25.*z-1.)*std::log(1.-z)/36.+std::log(1.-z)/z/2.-(1.+z)*PI2/36.+(-49.+38.*z*z-55.*z)/72.;
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

double BXsDecay::h88(double delta) {
    return 4./27.*(((1.+0.5*delta)*delta*log(delta)-6.*log(1.-delta)-2.*Li2(1.-delta)+PI2/3.-16./3.*delta-5./3.*delta*delta+delta*delta*delta/9.)*log(cache.m_b_1S/cache.m_s)-2.*Li3(delta)+(5.-2.*log(delta))*(Li2(1.-delta)-PI2/6.)-PI2/12.*delta*(2.+delta)+(0.5*delta+0.25*delta*delta-log(1.-delta))*pow(log(delta),2.)+(151./18.-PI2/3.)*log(1.-delta)+(-53./12.-19./12.*delta+2./9.*delta*delta)*delta*log(delta)+787./72.*delta+227./72.*delta*delta-41./72.*delta*delta*delta);
}

double BXsDecay::h77(double delta) {
    auto f = [this] (double z) {
        return F2nf(z);
    };

    return 4 * integrate(f, 0, 1 - delta, 1e-3);
}

std::array<std::array<double, 8>, 8> BXsDecay::phi_2_b0(double delta, double z) {
    std::array<std::array<double, 8>, 8> phi {};
    std::array<std::array<double, 8>, 8> phi_ij_1 = phi_1(delta, z);
    phi[1][1] = cache.beta_0 * (phi_ij_1[1][1] * cache.L_b + h22(z, delta));
    phi[1][6] = cache.beta_0 * (phi_ij_1[1][6] * cache.L_b + h27(z, delta));
    phi[1][7] = cache.beta_0 * (phi_ij_1[1][7] * cache.L_b + h28(z, delta));
    phi[6][6] = cache.beta_0 * (phi_ij_1[6][6] * cache.L_b + h77(z));
    phi[7][7] = cache.beta_0 * (phi_ij_1[7][7] * cache.L_b + h88(z));
    phi[0][0] = phi[1][1] / 36.0;
    phi[0][1] = -phi[1][1] / 3.0;
    phi[0][6] = -phi[1][6] / 6.0;
    phi[0][7] = -phi[1][7] / 6.0;
    return phi;
}

std::array<double, 8> BXsDecay::r_hat_2(double z) {
    std::array<double, 8> r {};
    double Lb = cache.L_b;
    double Lb2 = Lb * Lb;
    r[1] = -3. / 2 * r22(z) + 2 * (a(z) + b(z) - 290. / 81.) * Lb - 100. / 81 * Lb2;
    r[0] = -1. / 6 * r[1];
    r[6] = -3803. / 54 - 46 * PI2 / 27 + 80 * ZETA3 / 3 + (8 * PI2 / 9 - 98. / 3) * Lb - 16 * Lb2 / 3;
    r[7] = 1256. / 81 - 64 * PI2 / 81 - 32 * ZETA3 / 9 + (-8 * PI2 / 27 + 188. / 27) * Lb + 8 * Lb * Lb / 9;
    return r;
}

std::array<std::array<double, 8>, 8> BXsDecay::K_2_b0() {
    std::array<std::array<double, 8>, 8> K {};
    auto r = r_hat_2(cache.z);
    auto phi = phi_2_b0(cache.delta, cache.z);

    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 8; j++) {
            K[i][j] = K[j][i] 
                    = 2 * (1 + kron(i, j)) * phi[i][j] + cache.beta_0 * kron(j, 6) * (r[i]);
        }
    }

    return K;
}

double BXsDecay::F2a(double z) {
    if(fpeq(z, 0.)) return 0;

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

double BXsDecay::phi_77_rem(double delta) {
    auto f = [this] (double z) {
        return (16 * F2a(z) + 36 * F2na(z) + 87 * F2nf(z)) / 9.; 
    };

    double phi_77_A = integrate(f, 0, 1 - delta, 1e-3);
    double ld = std::log(delta);
    return -4 * phi_77_A - 8 * PI * cache.alpha_s_upsilon * (2 * delta * ld * ld + (4 + delta * (7 + delta * (-2 + delta))) * ld + 7 + delta * (-8. / 3 + delta * (-7 + delta * (4 - 4 * delta / 3)))) / (27 * delta);
}

std::array<std::array<double, 8>, 8> BXsDecay::K_2_rem(double z) {
    std::array<std::array<double, 8>, 8> K {};
    std::array<std::array<double, 8>, 8> K_1 = this->K_1();
    double L_D = cache.L_b - std::log(z);
    K[1][1] = std::pow((218. / 243 - 208. * L_D / 81), 2);
    K[0][0] = K[1][1] / 36.0;
    K[0][1] = K[1][0] = -K[1][1] / 6.0;
    K[1][6] = K[6][1] = (218./243.-208./81.*L_D) * K_1[6][6] + (127. / 324 - 35. * L_D / 27) * K_1[6][7] + 2 * (1 - L_D) * (K_1[3][6] - cache.beta_0 * (26. / 81 - 4 * cache.L_b / 27)) / 3 + L_D * (1150 - 4736 * L_D) / 729 - 1617980. / 19683 + 20060 * ZETA3 / 243 + 1664 * cache.L_c / 81;
    K[1][7] = K[7][1] = (218./243.-208./81.*L_D) * K_1[6][7] + (127. / 324 - 35. * L_D / 27) * K_1[7][7] + 2. * (1 - L_D) * K_1[3][7] / 3.;
    K[0][6] = K[6][0] = -K[1][6] / 6 + (5. / 16 - 3 * L_D / 4) * K_1[6][7] - 1237. / 729 + 232 * ZETA3 / 27 + L_D * (-20 + 70 * L_D) / 27;
    K[0][7] = K[7][0] = -K[1][7] / 6 + (5. / 16 - 3 * L_D / 4) * K_1[7][7];
    K[6][6] = (K_1[6][6] - 4 * phi_77(cache.delta) + 2 * std::log(z) / 3) * K_1[6][6] + L_D * (224. / 27 - 32 * L_D / 9) - 79.2838955662 + cache.L_b * (256 * PI2 / 27 - 2720. / 9 - 160 * cache.L_b / 3) + 512 * PI * cache.alpha_s_upsilon / 27 + 4 * phi_77_rem(cache.delta);
    K[6][7] = K[7][6] = (-50. / 3 + 8 * PI2 / 3 - 2 * L_D / 3) * K_1[6][7] + L_D * (-112. / 81 + 16 * L_D / 27) + 364. / 243;
    K[7][7] = (-50. / 3 + 8 * PI2 / 3 - 2 * L_D / 3) * K_1[7][7];

    return K;
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

double BXsDecay::r22_large_z(double z) {
    double lz = std::log(z);
    return 27650. / 6561 + lz * (112. / 243 + 8 * lz / 9);
}

double BXsDecay::P22_rem() {
    complex_t r21_0 {-1666. / 243, -80 * PI / 81};
    double r22_0 {67454. / 6561 - 124 * PI2 / 729}; 
    double x1 = std::pow(std::abs(cache.C_b_LO[WCoef::C1]), 2) / 36 + std::pow(std::abs(cache.C_b_LO[WCoef::C2]), 2) - std::real(cache.C_b_LO[WCoef::C1] * std::conj(cache.C_b_LO[WCoef::C2])) / 3
                + std::pow(std::abs(cache.C_b_LO[WCoef::CP1]), 2) / 36 + std::pow(std::abs(cache.C_b_LO[WCoef::CP2]), 2) - std::real(cache.C_b_LO[WCoef::CP1] * std::conj(cache.C_b_LO[WCoef::CP2])) / 3;
    double x2 = std::real(cache.C_b_LO[WCoef::C7] * std::conj(4019. * cache.C_b_LO[WCoef::C1] / 486. - 1184. * cache.C_b_LO[WCoef::C2] / 81. - 4. * cache.C_b_LO[WCoef::C7] + 4. * cache.C_b_LO[WCoef::C8] / 3.))
                + std::real(cache.C_b_LO[WCoef::CP7] * std::conj(4019. * cache.C_b_LO[WCoef::CP1] / 486. - 1184. * cache.C_b_LO[WCoef::CP2] / 81. - 4. * cache.C_b_LO[WCoef::CP7] + 4. * cache.C_b_LO[WCoef::CP8] / 3.));

    double phi_77_1 = phi_77(cache.delta);
    double K1_77 = 4 * phi_77_1 * (-182. / 9 + 8 * PI2 / 9 - gamma_i7[6] * cache.L_b);
    double K77rem_z0 = (K1_77-4.*phi_77_1+2./3.*cache.L_b)*K1_77-587708./729.-628./405.*pow(PI,4.)
	+32651./729.*PI2+428./27.*PI2*log(2.)+25150./81.*ZETA3-448./9.*cache.L_b*cache.L_b+(80./9.*PI2-2524./9.)*cache.L_b
	+512./27.*PI*cache.alpha_s_upsilon+4.*phi_77_rem(cache.delta)-8.*(phi_77_1 * cache.L_b + h77(cache.z))/3.;
    double x5 = K77rem_z0 * (std::pow(std::abs(cache.C_b_LO[WCoef::C7]), 2) + std::pow(std::abs(cache.C_b_LO[WCoef::CP7]), 2));

    auto target = [this, r21_0, r22_0, x1, x2, x5] (double z) {
        auto K_2 = K_2_rem(z);
        double y = gen_P00(K_2);
        // printf("P22rem = %.5e\n", y);
        double a1 = std::pow(std::abs(r2_large_z(z)), 2) - std::pow(std::abs(r21_0), 2);
        double a2 = r22_large_z(z) - r22_0;
        return y - a1 * x1 - a2 * x2 - x5;
    };

    double a_3_z0 = r2_large_z(cache.z0) - std::real(r21_0); 
    double a_3_z1 = r2_large_z(cache.z1) - std::real(r21_0);
    double a_4 = 2.*(4./3.-4./81.);
    double y_0 = target(cache.z0);
    double y_1 = target(cache.z1);

    // printf("a3 = %.5e\n", a_3_z0);
    // printf("b3 = %.5e\n", a_3_z1);
    // printf("a4 = %.5e\n", a_4);
    // printf("b4 = %.5e\n", a_4);
    // printf("a = %.5e\n", y_0);
    // printf("b = %.5e\n", y_1);

    double x3 = (y_1 - y_0) / (a_3_z1 - a_3_z0);
    double x4 = (y_0 * a_3_z1 - y_1 * a_3_z0) / (a_4 * (a_3_z1 - a_3_z0));
    complex_t r_21 = -1666. / 243 + 2 * (a(cache.z) + b(cache.z)) - 80. * I * PI / 81.;

    // printf("x1 = %.5e\n", x1);
    // printf("x2 = %.5e\n", x2);
    // printf("x3 = %.5e\n", x3);
    // printf("x4 = %.5e\n", x4);
    // printf("x5 = %.5e\n", x5);

    return x1 * (std::pow(std::abs(r_21), 2) - std::pow(std::abs(r21_0), 2))
            + x2 * (r22(cache.z) - r22_0)
            + x3 * std::real(r_21 - r21_0)
            + x4 * dr2_dlogz(cache.z)
            + x5;
}

double BXsDecay::P() {
    double p0 = std::pow(std::abs(cache.C_b_LO[WCoef::C7]), 2) + std::pow(std::abs(cache.C_b_LO[WCoef::CP7]), 2);
    double p11 = 2 * std::real(cache.C_b_LO[WCoef::C7] * std::conj(cache.C_b_NLO[WCoef::C7]));
    double p12 = std::pow(std::abs(cache.C_b_NLO[WCoef::C7]), 2) + 2 * std::real(cache.C_b_LO[WCoef::C7] * std::conj(cache.C_b_NNLO[WCoef::C7]));
    double p21 = gen_P00(K_1());
    double p32 = gen_P01(K_1());
    double p22 = gen_P00(K_2_b0()) + P22_rem();

    // LOG_INFO("p0", p0);
    // LOG_INFO("p11", p11);
    // LOG_INFO("p12", p12);
    // LOG_INFO("p21", p21);
    // LOG_INFO("p32", p32);
    // LOG_INFO("p22", p22);
    
    double k = cache.alpha_s_mu_b / (4 * PI);
    return p0 + k * ((p11 + p21) + k * (p12 + p22 + p32));
}

double BXsDecay::N() {
    double Kc = 0.0;
    for (size_t i = 0; i < a_i.size(); i++)
        Kc += d_i[i] * std::pow(cache.eta, a_i[i]);

    double Kt = std::real(
        (cache.C_w[WCoef::C7] + 23. / 36) * std::pow(cache.eta, 4.0 / 23.0)
        - 8. * (cache.C_w[WCoef::C8] + 1. / 3) * (std::pow(cache.eta, 4.0 / 23.0) - std::pow(cache.eta, 2.0 / 23.0)) / 3.
    );
    double eta_factor = std::pow(cache.eta, 6.0 / 23.0) + std::pow(cache.eta, -12.0 / 23.0);
    return -(Kc + cache.r_msbar_1S * Kt) * eta_factor * cache.lambda_2 / (18 * std::pow(cache.m_c, 2));
}

double BXsDecay::C2_em(double eta) {
    return -190 * std::pow(eta, -35. / 23) / 8073 - 359 * std::pow(eta, -17. / 23) / 3105 + 4276 * std::pow(eta, -12. / 23) / 121095
            + 350531 * std::pow(eta, -9. / 23) / 1009125 + 2 * std::pow(eta, -7. / 23) / 4347 - 5956 * std::pow(eta, 6. / 23) / 15525
            + 38380 * std::pow(eta, 14. / 23) / 169533 - 748 * std::pow(eta, 16. / 23) / 8625;
}

double BXsDecay::C8_em(double eta) {
    return -32 * std::pow(eta, -9. / 23) / 575 + 32 * std::pow(eta, -7. / 23) / 1449 + 640 * std::pow(eta, 14. / 23) / 1449 - 704 * std::pow(eta, 16. / 23) / 1725;
}

scalar_t BXsDecay::C7_em(double eta) {
    return (32 * std::pow(eta, -9. / 23) / 75 - 40 * std::pow(eta, -7. / 23) / 69 + 88 * std::pow(eta, 16. / 23) / 575) * cache.C_w[WCoef::C7] + C8_em(eta) * cache.C_w[WCoef::C8] + C2_em(eta);
}

double BXsDecay::epsilon_em() {
    double k_SL = 2. * cache.alpha_s_mu_b * std::log(cache.mu_W / cache.mu_b) / PI;
    // TODO : check eta or alpha_mub
    return (2 * std::real(C7_em(cache.eta) * std::conj(cache.C_b_LO[WCoef::C7])) - k_SL * std::pow(std::abs(cache.C_b_LO[WCoef::C7]), 2)) * cache.alpha_em / cache.alpha_s_mu_b;
}

double BXsDecay::C() {
    double delta_as = QCDHelper::alpha_s(4.6, MassType::MSBAR) - 0.22;
    double delta_b = cache.m_b_mb - 4.18;
    double delta_c = cache.m_c_3gev - 1.0;
    double rho = std::pow(cache.m_c_3gev / cache.m_b_1S, 2);
    double g = 1. + rho * (-8. + rho * (-12. * std::log(rho) + rho * (8. - rho)));
    return g * (0.849 - 0.92 * delta_as + 0.0596 * delta_b - 0.2237 * delta_c - 0.0167 * cache.mu_G2 - 0.203 * cache.rho_D3 + 0.004 * cache.rho_LS3);
}

double BXsDecay::BR_B_Xs_gamma() {
    double p = this->P();
    double n = this->N();
    double epsilon_em = this->epsilon_em();
    
    return cache.BR_B__Xc_e_nu_exp * cache.ckm_factor * 6 * cache.alpha_em / (PI * C()) * (p + n + epsilon_em);
}

std::vector<ObservableValue> BXsDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_B_XS_GAMMA:   
        value = BR_B_Xs_gamma();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BXsDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}