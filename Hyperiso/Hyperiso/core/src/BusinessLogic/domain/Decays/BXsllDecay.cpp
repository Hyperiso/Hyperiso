#include "BXsllDecay.h"

void BXsllDecay::load_params() {
    cache.cc_res_mass = {3.096916, 3.68609,   3.77292, 4.039 , 4.153  , 4.421 };
    cache.cc_res_br = {5.93e-2 , 7.7e-3 ,   1.1e-5 , 1.4e-5, 1.0e-5 , 1.1e-5};
    cache.cc_res_width_tot = {9.29e-5 , 3.04e-4,   2.73e-2, 8.0e-2, 1.03e-1, 6.2e-2};
    cache.cc_res_width_had = {8.147e-5, 2.9746e-4, 2.36e-2, 5.2e-2, 7.8e-2 , 4.3e-2};
    cache.q2_low_bound = {1., 6.};
    cache.q2_high_bound = {14.4, 22.};

    ObsParameterProxy p;
    cache.alpha_em = 1. / p(ParamId{ParameterType::SM, "SMINPUTS", 1});
    double m_c = p(ParamId{ParameterType::SM, "MASS", 4});
    cache.m_b_1S = p(ParamId{ParameterType::SM, "QCD", {5, 3}});
    complex_t V_tb = p(ParamId{ParameterType::SM, "VCKM", {2, 2}});
    complex_t V_ts = p(ParamId{ParameterType::SM, "VCKM", {2, 1}});
    complex_t V_cb = p(ParamId{ParameterType::SM, "VCKM", {1, 2}});
    complex_t V_cs = p(ParamId{ParameterType::SM, "VCKM", {1, 1}});
    cache.m_mu_hat = p(ParamId{ParameterType::SM, "MASS", 13}) / cache.m_b_1S;
    cache.m_tau_hat = p(ParamId{ParameterType::SM, "MASS", 15}) / cache.m_b_1S;
    cache.m_c_hat = m_c / cache.m_b_1S;
    cache.m_D_hat = p(ParamId{ParameterType::FLAVOR, "FMASS", 421}) / cache.m_b_1S;
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(w_config.hadronic_scale, MassType::POLE, MassType::POLE));
    cache.z = std::pow(cache.m_c_hat, 2);
    cache.L_b = std::log(w_config.hadronic_scale / cache.m_b_1S);
    cache.L_b_5GeV = std::log(w_config.hadronic_scale / 5.0);
    cache.L_l_mu = -2 * std::log(cache.m_mu_hat);
    cache.L_l_tau = -2 * std::log(cache.m_tau_hat);

    cache.pref_dB0_ds = 1. + 3. * p(ParamId{ParameterType::DECAY, "B_Xs", 6}) * g_lambda(cache.z) / (2. * pow(cache.m_b_1S, 2) * f(cache.z)) - p(ParamId{ParameterType::DECAY, "B_Xsll", 2}) * g_rho(cache.z) / (6. * pow(cache.m_b_1S, 3) * f(cache.z));
    cache.pref_dB_ds = p(ParamId{ParameterType::DECAY, "B_Xs", 2}) * std::pow(std::abs(V_tb * std::conj(V_ts) / V_cb), 2) / (std::pow(2 * PI / cache.alpha_em, 2) * f(cache.z) * kappa(cache.z));
    cache.pref_A0_0 = 1 + 3 * p(ParamId{ParameterType::DECAY, "B_Xs", 6}) * g_lambda(cache.z) / (2 * std::pow(cache.m_b_1S, 2) * f(cache.z));
    cache.pref_A0_1 = 4 * p(ParamId{ParameterType::DECAY, "B_Xsll", 1}) / (3 * std::pow(cache.m_b_1S, 2));
    cache.pref_delta_mb2 = 3. * p(ParamId{ParameterType::DECAY, "B_Xs", 6}) / (2. * std::pow(cache.m_b_1S * std::abs(V_tb), 2));
    cache.pref_delta_mb3 = -p(ParamId{ParameterType::DECAY, "B_Xsll", 2}) / (pow(cache.m_b_1S, 3) * std::pow(std::abs(V_tb), 2));
    cache.pref_delta_mc2 = 8. * p(ParamId{ParameterType::DECAY, "B_Xs", 6}) / (9. * std::pow(m_c, 2)) * std::abs(std::conj(V_cs) * V_cb / (std::conj(V_ts) * std::pow(V_tb, 3)));
    cache.pref_delta_brems = cache.alpha_s_mu_b / (4. * PI);
    cache.pref_delta_em = cache.alpha_em / (4. * PI);

    cache.C = w_proxy->getAFR(WGroup::B, w_config.order);
    auto C_P = w_proxy->getAFR(WGroup::BPrime, w_config.order);
    auto C_S = w_proxy->getAFR(WGroup::BScalar, w_config.order);
    cache.C.insert(C_P.begin(), C_P.end());
    cache.C.insert(C_S.begin(), C_S.end());

    cache.C_LO = w_proxy->getAR(WGroup::B, QCDOrder::LO);
    auto C_P_LO = w_proxy->getAR(WGroup::BPrime, QCDOrder::LO);
    cache.C_LO.insert(C_P_LO.begin(), C_P_LO.end());

    size_t ff_order {20};
    fill_cache(f_17, 0, 1, cache.F_17_lookup, cache.L_b, cache.z, ff_order);
    fill_cache(f_27, 0, 1, cache.F_27_lookup, cache.L_b, cache.z, ff_order);
    fill_cache(f_19_1S, 0, 1, cache.F_19_lookup, cache.L_b, cache.z, ff_order);
    fill_cache(f_29_1S, 0, 1, cache.F_29_lookup, cache.L_b, cache.z, ff_order);

    auto hatify = [this] (double q2) { return q2 / pow(cache.m_b_1S, 2); };
    cache.s_hat_low_bound = {hatify(cache.q2_low_bound.first), hatify(cache.q2_low_bound.second)};
    cache.s_hat_high_bound = {hatify(cache.q2_high_bound.first), hatify(cache.q2_high_bound.second)};

    auto bound_func = std::bind(&BXsllDecay::delta_bremB_base, &*this, std::placeholders::_1);
    fill_cache(bound_func, cache.s_hat_low_bound.first, cache.s_hat_low_bound.second, cache.delta_brems_lookup_low);
    fill_cache(bound_func, cache.s_hat_high_bound.first, cache.s_hat_high_bound.second, cache.delta_brems_lookup_high);
}

double BXsllDecay::f(double z) {
    return 1 - 8 * z + 8 * pow(z, 3) - pow(z, 4) - 12 * pow(z, 2) * log(z);
}

double BXsllDecay::h(double z) {
    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;
    double lz = log(z);
    double lw = log(1 - z);
    double lz2 = lz * lz;

    return -(1 - z2) * (25./4 - 239./3 * z + 25./4 * z2) + z * lz * (20. + 90. * z - 4./3 * z2 + 17./3 * z3)
           + z2 * lz2 * (36. + z2) + (1. - z2) * (17./3. - 64./3. * z + 17./3. * z2) * lw - 4. * (1. + 30. * z2 + z4) * lz * lw
	       - (1. + 16. * z2 + z4) * (6. * Li2(z) - PI2) - 32. * pow(z, 1.5) * (1. + z) * (PI2 - 4. * Li2(sqrt(z)) + 4. * Li2(-sqrt(z)) - 2. * lz * log((1.-sqrt(z))/(1.+sqrt(z))));
}

double BXsllDecay::g_rho(double z) {
    return 77.-88.*z+24.*pow(z,2.)-8.*pow(z,3)+5.*pow(z,4)+48.*log(z)+36.*pow(z,2)*log(z);
}

double BXsllDecay::g_lambda(double z) {
    return 3.-8.*z+24.*pow(z,2)-24.*pow(z,3)+5.*pow(z,4)+12.*pow(z,2)*log(z);
}

double BXsllDecay::kappa(double z) {
    return 1 - 2. * cache.alpha_s_mu_b * h(z) / (3. * PI * f(z));
}

double BXsllDecay::f_7(double s) {
    return 1./6./(s-1.)/(s-1.)*(24.*(1.+13.*s-4.*s*s)*Li2(sqrt(s))+12.*(1.-17.*s+6.*s*s)*Li2(s)+6.*s*(6.-7.*s)*log(s)
	+24.*(1.-s)*(1.-s)*log(s)*log(1.-s)+12.*(-13.+16.*s-3.*s*s)*(log(1.-sqrt(s))-log(1.-s))
	+39.-2.*PI2+252.*s-26.*PI2*s+21.*s*s+8.*PI2*s*s-180.*sqrt(s)-132.*s*sqrt(s));
}

double BXsllDecay::f_9(double s) {
    return -1./6./(s-1.)/(s-1.)*(48.*s*(-5.+2.*s)*Li2(sqrt(s))+24.*(-1.+7.*s-3.*s*s)*Li2(s)+6.*s*(-6.+7.*s)*log(s)
	-24.*(1.-s)*(1.-s)*log(s)*log(1.-s)+24.*(5.-7.*s+2.*s*s)*(log(1.-sqrt(s))-log(1.-s))
	-21.-156.*s+20.*PI2*s+9.*s*s-8.*PI2*s*s+120.*sqrt(s)+48.*s*sqrt(s));
}

complex_t BXsllDecay::Gm1(double t) {
    complex_t L = log((sqrt((complex_t)t)+sqrt((complex_t)t-4.))/2.);
    return -2.*I*PI*L-PI2/2.+2.*pow(L,2.);
}

complex_t BXsllDecay::G0(double t) {
    complex_t rt = sqrt(((complex_t)t-4.)/(complex_t)t);
    return -I*PI*rt-2.+2.*rt*log((sqrt((complex_t)t)+sqrt((complex_t)t-4.))/2.);
}

complex_t BXsllDecay::Delta_i_23(double s, double z, double w) {
    return -2.+4./(w-s)*(z*Gm1(s/z)-z*Gm1(w/z)-s/2.*G0(s/z)+s/2.*G0(w/z));
}

complex_t BXsllDecay::Delta_i_27(double s, double z, double w) {
    return 2.*(G0(s/z) - G0(w/z));
}

double BXsllDecay::tau_22(double s, double w, complex_t Delta_23, complex_t Delta_27) {
    return 8./27.*(w-s)*(1.-w)*(1.-w)/s/w/w/w*((3.*w*w+2.*s*s*(2.+w)-s*w*(5.-2.*w))*pow(abs(Delta_23),2.)
	+(2.*s*s*(2.+w)+s*w*(1.+2.*w))*pow(abs(Delta_27),2.)
	+4.*s*(w*(1.-w)-s*(2.+w))*real(Delta_23*conj(Delta_27)));
}

complex_t BXsllDecay::tau_27(double s, double w, complex_t Delta_23, complex_t Delta_27) {
    return 8./3./s/w*(((1.-w)*(4.*s*s-s*w+w*w)+s*w*(4.+s-w)*log(w))*Delta_23
	-4.*s*s*(1.-w)+s*w*(4.+s-w)*log(w)*Delta_27);
}

complex_t BXsllDecay::tau_28(double s, double w, complex_t Delta_23, complex_t Delta_27) {
    return 8./9./s/w/(w-s)*((pow(w-s,2.)*(2.*s-w)*(1.-w))*Delta_23
	-(2.*s*pow(w-s,2.)*(1.-w))*Delta_27
	+s*w*((1.+2.*s-2.*w)*Delta_23-2.*(1.+s-w)*Delta_27)*log(s/((1.+s-w)*(w*w+s*(1.-w)))));
}

complex_t BXsllDecay::tau_29(double s, double w, complex_t Delta_23, complex_t Delta_27) {
    return 4./3./w*((2.*s*(1.-w)*(s+w)+4.*s*w*log(w))*Delta_23
	-(2.*s*(1.-w)*(s+w)+w*(3.*s+w)*log(w))*Delta_27);
}

double BXsllDecay::tau_77(double s) {
    return -2./9./(2.+s)*(2.*(1.-s)*(1.-s)*log(1.-s)+6.*s*(2.-2.*s-s*s)/(1.-s)/(1.-s)*log(s)+(11.-7.*s-10.*s*s)/(1.-s));
}

double BXsllDecay::tau_78(double s) {
    return 8./9./s*(25.-2.*PI2-27.*s+3.*s*s-s*s*s+12.*(s+s*s)*log(s)+6.*pow((PI/2.-atan((2.-4.*s+s*s)/(2.-s)*sqrt(s)*sqrt(4.-s))),2.)
	-24.*real(CLi2((s-I*sqrt(s)*sqrt(4.-s))/2.))-12.*((1.-s)*sqrt(s)*sqrt(4.-s)-atan((sqrt(s)*sqrt(4.-s))/(2.-s)))
	*(atan(sqrt((4.-s)/s))-atan((sqrt(s)*sqrt(4.-s))/(2.-s))));
}

double BXsllDecay::tau_88(double s) {
    return 4./27./s*(-8.*PI2+(1.-s)*(77.-s-4.*s*s)-24.*Li2(1.-s)+3.*(10.-4.*s-9.*s*s+8.*log((sqrt(s))/(1.-s)))*log(s)
	+48.*real(CLi2((3.-s)/2.+I*(1.-s)*sqrt(4.-s)/2./sqrt(s)))-6.*((20.*s+10.*s*s-3.*s*s*s)/sqrt(s)/sqrt(4.-s)-8.*PI+8.*atan(sqrt((4.-s)/s)))
	*(atan(sqrt((4.-s)/s))-atan((sqrt(s)*sqrt(4.-s))/(2.-s))));
}

double BXsllDecay::tau_89(double s) {
    return 2./3.*(s*(4.-s)-3.-4.*log(s)*(1.-s-s*s)
	-8.*real(CLi2(s/2.+I*sqrt(s)*sqrt(4.-s)/2.)-CLi2((-2.+s*(4.-s))/2.+I*((2.-s)*sqrt(s)*sqrt(4.-s))/2.))
	+4.*(s*s*sqrt((4.-s)/s)+2.*atan(sqrt(s)*sqrt(4.-s)/(2.-s)))*(atan(sqrt((4.-s)/s))-atan((sqrt(s)*sqrt(4.-s))/(2.-s))));
}

double BXsllDecay::tau_99(double s) {
    return -4./9./(1.+2.*s)*(2.*(1.-s)*(1.-s)*log(1.-s)+3.*s*(1.+s)*(1.-2.*s)/(1.-s)/(1.-s)*log(s)+3.*(1.-3.*s*s)/(1.-s));
}

double BXsllDecay::tau_79(double s) {
    return -4.*(1.-s)*(1.-s)/9./s*log(1.-s)-4.*s*(3.-2.*s)*log(s)/9./(1.-s)/(1.-s)-2./9.*(5.-3.*s)/(1.-s);
}

complex_t BXsllDecay::tau_210(double s, double z) {
    auto f = [&] (double w) {
        return -s/(s-w)/(1.-s)/(1.-s)*((4.*(1.-s)*(1.+w)-2.*fabs(s-w*w)*(w*(3.+w)-s*(1.-w))/w/w
            +(2.+5.*w+2.*w*w+s*(3.+4.*w))*log((s+w*w+fabs(s-w*w))/2./w)-(s-w)/s/sqrt((1.+w)*(1.+w)-4.*s)
            *(w*(2.-w)-s*(6.-5.*w))*(log(1.+w-s*(3.-w)+(1.-s)*sqrt((1.+w)*(1.+w)-4.*s))
            -log(s*(1.-3.*w)+w*w*(1.+w)+fabs(s-w*w)*sqrt((1.+w)*(1.+w)-4.*s))))*Delta_i_23(s, z, w)
            -(2.*(1.-s)*(1.+2.*w)-2.*fabs(s-w*w)*(w*(2.+w)-s*(1.-w))/w/w
            +2.*(s*(1.+2.*w)+w*(2.+w))*log((s+w*w+fabs(s-w*w))/2./w)
            +4.*(1.-w)*(s-w)/sqrt((1.+w)*(1.+w)-4.*s)*(log(1.+w-s*(3.-w)+(1.-s)*sqrt((1.+w)*(1.+w)-4.*s))
            -log(s*(1.-3.*w)+w*w*(1.+w)+fabs(s-w*w)*sqrt((1.+w)*(1.+w)-4.*s))))*Delta_i_27(s, z, w));
    };

    return c_integrate(f, s, 1, 1e-2);
}

double BXsllDecay::tau_710(double s) {
    return -5./2.+1./3./(1.-3.*s)-1./3.*s*(6.-7.*s)*log(s)/(1.-s)/(1.-s)-1./9.*(3.-7.*s+4.*s*s)*log(1.-s)/s+f_7(s)/3.;
}

double BXsllDecay::tau_810(double s){
    return 1./6./(1.-s)/(1.-s)*(3.*((1.-sqrt(s))*(1.-sqrt(s))*(23.-6.*sqrt(s)-s)+4.*(1.-s)*(7.+s)*log(1.+sqrt(s))
	+2.*s*(1.+s-log(s))*log(s))+2.*(-3.*PI2*(1.+2.*s)+6.*(3.-s)*s*log(2.-sqrt(s))
	-36.*(1.+2.*s)*CLi2(-sqrt(s))-6.*sqrt(s/(4.-s))*(2.*(-3.+s)*s*atan((2.+sqrt(s))/sqrt(4.-s))
	+2.*PI*log(2.-sqrt(s))-atan(sqrt((4.-s)/s))*((-3.+s)*s+4.*log(2.-sqrt(s)))
	-atan(sqrt(s*(4.-s))/(2.-s))*((-3.+s)*s-log(s))+4.*real(I*CLi2((-2.+I*sqrt(4.-s)+sqrt(s))*sqrt(s)/(I*sqrt(4.-s)-sqrt(s))))
	-2.*real(I*CLi2(I/2.*sqrt(4.-s)*(1.-s)*sqrt(s)+(3.-s)*s/2.)))));
}

double BXsllDecay::tau_910(double s) {
    return -5./2.+1./3./(1.-s)-1./3.*s*(6.-7.*s)*log(s)/(1.-s)/(1.-s)-2./9.*(3.-5.*s+2.*s*s)*log(1.-s)/s+f_9(s)/3.;
}

double BXsllDecay::sigma(double s) {
    return -4./3.*Li2(s)-2./3.*log(s)*log(1.-s)-2./9.*PI2-log(1.-s)-2./9.*(1.-s)*log(1.-s);
}

double BXsllDecay::sigma_9(double s) {
    return sigma(s)+1.5;
}

double BXsllDecay::sigma_7(double s, double L_mu) {
    return sigma(s)+1./6.-8./3.*L_mu;
}

complex_t BXsllDecay::F(double r) {
    if(r < 1) {
        return 3./2./r*(1./sqrt(r*(1.-r))*atan(sqrt(r/(1.-r)))-1.);
    } 
    return 3./2./r*(1./2./sqrt(r*(r-1.))*(log((1.-sqrt(1.-1./r))/(1.+sqrt(1.-1./r)))+I*PI)-1.);
}

complex_t BXsllDecay::Sigma_1(double s) {
    if (abs(s) < 0.4)
        return {
            23.787-120.948*s+365.373*s*s-584.206*s*s*s,
            1.653+6.009*s-17.080*s*s+115.880*s*s*s
        };

    double d = 1 - s;
    return {
        -148.061*d*d+492.539*d*d*d-1163.847*pow(d,4.)+1189.528*pow(d,5.),
        -261.287*d*d+1170.856*d*d*d-2546.948*pow(d,4.)+2540.023*pow(d,5.)
    };
}

double BXsllDecay::Sigma_2(double s) {
    if (abs(s) < 0.4)
        return 11.488-36.987*s+255.330*s*s-812.388*s*s*s+1011.791*s*s*s*s;

    double d = 1 - s;
    return -221.904*d*d+900.822*d*d*d-2031.620*pow(d,4.)+1984.303*pow(d,5.);
}

complex_t BXsllDecay::Sigma_3(double s) {
    if (abs(s) < 0.4)
        return {
            109.311-846.039*s+2890.115*s*s-4179.072*s*s*s,
            4.606+17.650*s-53.244*s*s+348.069*s*s*s
        };

    double d = 1 - s;
    return {
        -298.730*d*d+828.0675*d*d*d-2217.6355*pow(d,4.)+2241.792*pow(d,5.),
        -528.759*d*d+2095.723*d*d*d-4681.843*pow(d,4.)+5036.677*pow(d,5.)
    };
}

complex_t BXsllDecay::Sigma_7(double s, double z) {
    if (abs(s) < 0.4) {
        double a = pow(4 * z, 2); 
        return {
            -0.259023-28.424*s+205.533*s*s-603.219*s*s*s+722.031*s*s*s*s,
            (-12.20658-215.8208*(s-a)+412.1207*(s-a)*(s-a))*(s-a)*(s-a)*(s>a)
        };
    }
        
    double d = 1 - s;
    return {
        77.0256*d*d-264.705*d*d*d+595.814*pow(d,4.)-610.1637*pow(d,5.),
        135.858*d*d-618.990*d*d*d+1325.040*pow(d,4.)-1277.170*pow(d,5.)
    };
}

double BXsllDecay::omega_22(double s, double L_l) {
    return L_l * (Sigma_2(s)/8./(1.-s)/(1.-s)/(1.+2.*s)+real(Sigma_1(s))/9./(1.-s)/(1.-s)/(1.+2.*s)*cache.L_b_5GeV)
	+64./81.*omega_1010(s, L_l) * cache.L_b_5GeV * cache.L_b_5GeV;
}

complex_t BXsllDecay::omega_27(double s, double L_l) {
    return L_l*(Sigma_3(s)/96./(1.-s)/(1.-s))+8./9.*omega_79(s, L_l)*cache.L_b_5GeV;
}

complex_t BXsllDecay::omega_29(double s, double L_l) {
    return L_l * (Sigma_1(s)/8./(1.-s)/(1.-s)/(1.+2.*s))+16./9.*omega_1010(s,L_l)*cache.L_b_5GeV;
}

complex_t BXsllDecay::omega_210(double s, double L_l, double z) {
    return L_l*(-Sigma_7(s, z)/24./s/(1.-s)/(1.-s))+8./9.*omega_910(s,L_l)*cache.L_b_5GeV;
}

double BXsllDecay::omega_77(double s, double L_l) {
    return L_l*(s/2./(1.-s)/(2.+s)+log(1.-s)-s*(-3.+2.*s*s)/2./(1.-s)/(1.-s)/(2.+s)*log(s));
}

double BXsllDecay::omega_79(double s, double L_l) {
    return L_l*(-0.5/(1.-s)+log(1.-s)+(-1.+2.*s-2.*s*s)/2./(1.-s)/(1.-s)*log(s));
}

double BXsllDecay::omega_710(double s, double L_l) {
    return L_l*((7.-16.*sqrt(s)+9.*s)/4./(1.-s)+log(1.-sqrt(s))+(1.+3.*s)/(1.-s)*log((1.+sqrt(s))/2.)-s*log(s)/(1.-s));
}

double BXsllDecay::omega_99(double s, double L_l) {
    return L_l*(-(1.+4.*s-8.*s*s)/6./(1.-s)/(1.+2.*s)+log(1.-s)-(1.-6.*s*s+4.*s*s*s)*log(s)/2./(1.-s)/(1.-s)/(1.+2.*s))
	-Li2(s)/9.+4./27.*PI2-(37.-3.*s-6.*s*s)/72./(1.-s)/(1.+2.*s)-((41.+76.*s)*log(1.-s))/36./(1.+2.*s)
	+((6.-10.*s-17.*s*s+14.*s*s*s)/18./(1.-s)/(1.-s)/(1.+2.*s)+17.*log(1.-s)/18.)*log(s)-(1.-6.*s*s+4.*s*s*s)*log(s)*log(s)/2./(1.-s)/(1.-s)/(1.+2.*s);
}

double BXsllDecay::omega_910(double s, double L_l) {
    return L_l*(-(5.-16.*sqrt(s)+11.*s)/4./(1.-s)+log(1.-sqrt(s))+(1.-5.*s)/(1.-s)*log((1.+sqrt(s))/2.)-(1.-3.*s)*log(s)/(1.-s));
}

double BXsllDecay::omega_1010(double s, double L_l) {
    return L_l*(-(1.+4.*s-8.*s*s)/6./(1.-s)/(1.+2.*s)+log(1.-s)-(1.-6.*s*s+4.*s*s*s)*log(s)/2./(1.-s)/(1.-s)/(1.+2.*s));
}

complex_t BXsllDecay::g(double z, double s) {
    double z2=z*z;

	if(s==0.) 
        return -4./9.*log(z2)+8./27.-4./9.;
	
	if(z==0.) 
        return 8./27.-4./9.*(log(s)-I*PI);
	
	if(4.*z2<s) {
        return -4./9.*log(z2)+8./27.+16./9.*z2/s -2./9.*sqrt(1.-4.*z2/s)*(2.+4.*z2/s)*(log((1.+sqrt(1.-4.*z2/s))/(1.-sqrt(1.-4.*z2/s)))-I*PI);
    } else if(4.*z2>s) {
        return -4./9.*log(z2)+8./27.+16./9.*z2/s -4./9.*sqrt(4.*z2/s-1.)*(2.+4.*z2/s)*atan(1./sqrt(4.*z2/s-1.));
    }
        
	else 
        return -4./9.*log(z2)+8./27.+16./9.*z2/s;
}

double BXsllDecay::breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had) {
    double m_V_hat = m_V / cache.m_b_1S;
    double gamma_tot_hat = gamma_tot / cache.m_b_1S;
    double gamma_had_hat = gamma_had / cache.m_b_1S;
    return br * gamma_tot_hat * gamma_had_hat / (pow(s - pow(m_V_hat, 2.), 2.) + pow(m_V_hat * gamma_tot_hat, 2.));
}

double BXsllDecay::R_cc_cont(double s) {
    return s > 0.6 ? (s > 0.69 ? 1.02 : 11.33 * s - 6.8) : 0;
}

double BXsllDecay::R_cc(double s) {
    double R_cc_res;

    for (size_t k = 0; k < cache.cc_res_mass.size(); k++) {
        R_cc_res += breit_wigner(s, cache.cc_res_mass[k], cache.cc_res_br[k], cache.cc_res_width_tot[k], cache.cc_res_width_had[k]);
    }

    return 9. * s / pow(cache.alpha_em, 2) * R_cc_res + R_cc_cont(s);
}

double BXsllDecay::PV_breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had) {
    double s_c = 4 * pow(cache.m_D_hat, 2);
    double m_V_hat = m_V / cache.m_b_1S;
    double m_V_hat2 = pow(m_V_hat, 2);
    double gamma_tot_hat = gamma_tot / cache.m_b_1S;
    double gamma_had_hat = gamma_had / cache.m_b_1S; 
    double den = (pow(s - pow(m_V_hat, 2.), 2.) + pow(m_V_hat * gamma_tot_hat, 2.));
    double B = breit_wigner(s, m_V, br, gamma_tot, gamma_had);
    return 9 * s / pow(cache.alpha_em, 2) * B * (0.5 * log(den / pow(s_c - s, 2)) + (s - m_V_hat2) / gamma_tot_hat * m_V_hat * (atan((s_c - m_V_hat2) / gamma_tot_hat * m_V_hat) - PI / 2));
}

double BXsllDecay::PV_R_cc_cont(double s) {
    double s_c = 4 * pow(cache.m_D_hat, 2);
    return 1 / 3. * ((11.33 * s - 6.8) * log(abs((0.69 - s) / (s_c - s))) - 1.02 * log(abs(0.69 - s)) - 6.8 * log(s_c) - 2.902017); 
}

double BXsllDecay::PV_R_cc(double s) {
    double PV_res;

    for (size_t k = 0; k < cache.cc_res_mass.size(); k++) {
        PV_res += PV_breit_wigner(s, cache.cc_res_mass[k], cache.cc_res_br[k], cache.cc_res_width_tot[k], cache.cc_res_width_had[k]);
    }

    return PV_res + PV_R_cc_cont(s);
}

complex_t BXsllDecay::g_ld(double z, double s) {
    // LOG_INFO(R_cc(s));
    return g(z, 0) + (s * PV_R_cc(s) + I * PI * R_cc(s)) / 3.;
}

complex_t BXsllDecay::C9_eff(double s, QCDOrder order, bool prime) {
    s = std::clamp(s, 1e-6, 1. - 1e-6);
    auto C = order == QCDOrder::LO ? cache.C_LO : cache.C;
    auto C_ids = WCoefMapper::get_group(prime ? WGroup::BPrime : WGroup::B);
    complex_t g_0 = g(0, s);
    complex_t g_1 = g(1, s);
    complex_t g_mc = g_ld(cache.m_c_hat, s);

    LOG_INFO("g_mc =", g_mc);

    return C[C_ids[8]]
	        +(-32./27.*C[C_ids[0]]-8./9.*C[C_ids[1]]-16./9.*C[C_ids[2]]+32./27.*C[C_ids[3]]-112./9.*C[C_ids[4]]+512./27.*C[C_ids[5]]) * cache.L_b
	        +4./3.*C[C_ids[2]]+64./9.*C[C_ids[4]]+64./27.*C[C_ids[5]]
	        +g_mc*(4./3.*C[C_ids[0]]+C[C_ids[1]]+6.*C[C_ids[2]]+60.*C[C_ids[4]])
	        +g_1*(-7./2.*C[C_ids[2]]-2./3.*C[C_ids[3]]-38.*C[C_ids[4]]-32./3.*C[C_ids[5]])
	        +g_0*(-1./2.*C[C_ids[2]]-2./3.*C[C_ids[3]]-8.*C[C_ids[4]]-32./3.*C[C_ids[5]]);
}

complex_t BXsllDecay::F_17(double s) {
    return lerp(s, cache.F_17_lookup);
}

complex_t BXsllDecay::F_27(double s) {
    return lerp(s, cache.F_27_lookup);
}

complex_t BXsllDecay::F_19(double s) {
    return lerp(s, cache.F_19_lookup);
}

complex_t BXsllDecay::F_29(double s) {
    return lerp(s, cache.F_29_lookup);
}

complex_t BXsllDecay::C7_new(double s, bool prime) {
    auto C_0 = cache.C_LO;
    auto C7_eff = cache.C[prime ? WCoef::CP7 : WCoef::C7];
    // LOG_INFO("s_hat =", s);
    // LOG_INFO("C7eff =", C7_eff);
    // LOG_INFO("sigma_7 =", sigma_7(s, cache.L_b));
    // LOG_INFO("F_17 =", F_17(s));
    // LOG_INFO("F_27 =", F_27(s));
    // LOG_INFO("F_87 =", f_87(s, cache.L_b));
    return (1.+cache.alpha_s_mu_b/PI*sigma_7(s, cache.L_b))*C7_eff
	-cache.alpha_s_mu_b/4./PI*(C_0.at(prime ? WCoef::CP1 : WCoef::C1)*F_17(s)+C_0.at(prime ? WCoef::CP2 :WCoef::C2)*F_27(s)+C_0.at(prime ? WCoef::CP8 :WCoef::C8)*f_87(s, cache.L_b));
}

complex_t BXsllDecay::C9_new(double s, bool prime) {
    if (abs(s - 1) < 1e-6) s = 1;
    auto C_0 = cache.C_LO;
    LOG_INFO("C9_eff =", C9_eff(s, this->w_config.order, prime));
    return (1.+cache.alpha_s_mu_b/PI*sigma_9(s))*C9_eff(s, this->w_config.order, prime)
	        -cache.alpha_s_mu_b/4./PI*(C_0.at(prime ? WCoef::CP1 :WCoef::C1)*F_19(s)+C_0.at(prime ? WCoef::CP2 :WCoef::C2)*F_29(s)+C_0.at(prime ? WCoef::CP8 :WCoef::C8)*f_89(s));
}

complex_t BXsllDecay::C10_new(double s, bool prime) {
    s = std::clamp(s, 1e-6, 1. - 1e-6);
    return (1.+cache.alpha_s_mu_b/PI*sigma_9(s))*(cache.C[prime ? WCoef::CP10 : WCoef::C10]);
}

double BXsllDecay::W_7(double s) {
    return pow(abs(C7_new(s, false)), 2) + pow(abs(C7_new(s, true)), 2);
}

double BXsllDecay::W_9(double s) {
    return pow(abs(C9_new(s, false)), 2) + pow(abs(C9_new(s, true)), 2);
}

double BXsllDecay::W_10(double s) {
    return pow(abs(C10_new(s, false)), 2) + pow(abs(C10_new(s, true)), 2);
}

complex_t BXsllDecay::W_27(double s) {
    return cache.C[WCoef::C2] * conj(C7_new(s, false)) + cache.C[WCoef::CP2] * conj(C7_new(s, true));
}

complex_t BXsllDecay::W_29(double s) {
    return cache.C[WCoef::C2] * conj(C9_new(s, false)) + cache.C[WCoef::CP2] * conj(C9_new(s, true));
}

complex_t BXsllDecay::W_210(double s) {
    return cache.C[WCoef::C2] * conj(C10_new(s, false)) + cache.C[WCoef::CP2] * conj(C10_new(s, true));
}

double BXsllDecay::W_79(double s) {
    return real(C7_new(s, false) * conj(C9_new(s, false))) + real(C7_new(s, true) * conj(C9_new(s, true)));
}

double BXsllDecay::W_710(double s) {
    return real(C7_new(s, false) * conj(C10_new(s, false))) + real(C7_new(s, true) * conj(C10_new(s, true)));
}

double BXsllDecay::W_910(double s) {
    return real(C9_new(s, false) * conj(C10_new(s, false)))  + real(C9_new(s, true) * conj(C10_new(s, true)));
}

double BXsllDecay::dB0_ds(double s, double ml_hat) {
    double W_Q1 = pow(abs(cache.C[WCoef::CQ1]), 2) + pow(abs(cache.C[WCoef::CPQ1]), 2);
    double W_Q2 = pow(abs(cache.C[WCoef::CQ2]), 2) + pow(abs(cache.C[WCoef::CPQ2]), 2);
    double W_10Q2 = real(cache.C[WCoef::CQ2] * conj(C10_new(s, false))) + real(cache.C[WCoef::CPQ2] * conj(C10_new(s, true)));

    double H7 = 4 * (1. + 2. * ml_hat * ml_hat / s) * (1. + 2. / s) * (1. + cache.alpha_s_mu_b / PI * tau_77(s));
    double H9 = (1. + 2. * ml_hat * ml_hat / s) * (1. + 2. * s) * (1. + cache.alpha_s_mu_b / PI * tau_99(s));
    double H10 = ((1.+2.*s)+2.*ml_hat*ml_hat/s*(1.-4.*s))*(1.+cache.alpha_s_mu_b/PI*tau_99(s));
    double H79 = 12 * (1.+2.*ml_hat*ml_hat/s)*(1.+cache.alpha_s_mu_b/PI*tau_79(s));

    LOG_INFO("W_7 =", W_7(s));
    LOG_INFO("W_9 =", W_9(s));
    LOG_INFO("W_10 =", W_10(s));
    LOG_INFO("W_79 =", W_79(s));

    return pow(1 - s, 2) * sqrt(1 - 4 * pow(ml_hat, 2) / s) * (
            H7 * W_7(s) + 
            H9 * W_9(s) + 
            H10 * W_10(s) + 
            H79 * W_79(s) +
            1.5 * (s - 4 * ml_hat * ml_hat) * W_Q1 +
            1.5 * s * W_Q2 +
            6 * ml_hat * W_10Q2
        );
}

double BXsllDecay::A_FB_0(double s, double ml_hat) {
    double W_7Q1 = real(C7_new(s, false) * cache.C[WCoef::CQ1]) + real(C7_new(s, true) * cache.C[WCoef::CPQ1]);
    double W_9Q1 = real(C9_new(s, false) * cache.C[WCoef::CQ1]) + real(C9_new(s, true) * cache.C[WCoef::CPQ1]);

    return pow(1 - s, 2) * sqrt(1 - 4 * pow(ml_hat, 2) / s) * (
            2 * (1 + cache.alpha_s_mu_b * tau_710(s) / PI) * W_710(s) + 
            s * (1 + cache.alpha_s_mu_b * tau_910(s) / PI) * W_910(s) + 
            ml_hat * (2 * W_7Q1 + W_9Q1)
        );
}

double BXsllDecay::delta_A_mb2(double s) {
    return s * (9 + 14 * s - 15 * pow(s, 2)) * W_910(s) + 2 * (7 + 10 * s - 9 * s * s) * W_710(s);
}

double BXsllDecay::delta_mb2(double s) {
    return -4 * (6 + 3 * s - 5 * pow(s, 3)) * W_7(s) / s + (1 - 15 * s * s + 10 * pow(s, 3)) * (W_9(s) + W_10(s)) - 4 * (5 + 6 * s - 7 * s * s) * W_79(s);
}

double BXsllDecay::delta_mb3(double s) {
    if (s > 0.4)
        return 0;

    return (5.*pow(s,4.)+19.*pow(s,3.)+9.*s*s-7.*s+22.)/6./(1.-s)*4.*W_7(s)/s+
            (10.*pow(s,4.)+23.*pow(s,3.)-9.*s*s+13.*s+11.)/6./(1.-s)*(W_9(s) + W_10(s))+
            4.*(-3.*pow(s,3.)+17.*s*s-s+3.)/2./(1.-s)* W_79(s);
}

double BXsllDecay::delta_mc2(double s) {
    complex_t f = F(s / (4. * cache.z));
    return pow(1 - s, 2) * real((1 + 6 * s - s * s) * f * W_27(s) / s + (2 + s) * f * W_29(s));
}

double BXsllDecay::delta_A_mc2(double s) {
    return pow(1 - s, 2) * real((1 + 3 * s) * F(s / (4. * cache.z)) * W_210(s));
}

double BXsllDecay::delta_bremA(double s) {
    complex_t C7_0 = cache.C_LO[WCoef::C7];
    complex_t C8_0 = cache.C_LO[WCoef::C8];
    complex_t CP7_0 = cache.C_LO[WCoef::CP7];
    complex_t CP8_0 = cache.C_LO[WCoef::CP8];
    complex_t C9_0 = C9_eff(s, QCDOrder::LO, false);
    complex_t CP9_0 = C9_eff(s, QCDOrder::LO, true);

    complex_t c_78 = ObsQCDProxy().get_constants()->C_F * (C7_0 * conj(C8_0) + CP7_0 * conj(CP8_0));
    complex_t c_88 = ObsQCDProxy().get_constants()->C_F * (C8_0 * conj(C8_0) + CP8_0 * conj(CP8_0));
    complex_t c_89 = ObsQCDProxy().get_constants()->C_F * (C8_0 * conj(C9_0) + CP8_0 * conj(CP9_0));

    return 2 * real(c_78 * tau_78(s) + c_89 * tau_89(s) + c_88 * tau_88(s));
}

double BXsllDecay::delta_bremB_base(double s) {
    auto C_0 = cache.C_LO;
    auto CP_0 = cache.C_LO;
    complex_t C9_0 = C9_eff(s, QCDOrder::LO, false);
    complex_t CP9_0 = C9_eff(s, QCDOrder::LO, true);

    double C_f = ObsQCDProxy().get_constants()->C_F;
    double C_tau_1 = C_f / (4 * pow(ObsQCDProxy().get_constants()->Nc, 2));
    double C_tau_2 = -C_f / (2 * ObsQCDProxy().get_constants()->Nc);

    complex_t c_11 = C_tau_1 * (C_0[WCoef::C1] * conj(C_0[WCoef::C1]) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP1]));
	complex_t c_12 = 2 * C_tau_2 * real((C_0[WCoef::C1]*conj(C_0[WCoef::C2])) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP2]));
	complex_t c_22 = C_f*(C_0[WCoef::C2]*conj(C_0[WCoef::C2]) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP2]));
	complex_t c_17 = C_tau_2*(C_0[WCoef::C1]*conj(C_0[WCoef::C7]) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP7]));
	complex_t c_27 = C_f*(C_0[WCoef::C2]*conj(C_0[WCoef::C7]) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP7]));
	complex_t c_18 = C_tau_2*(C_0[WCoef::C1]*conj(C_0[WCoef::C8]) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP8]));
	complex_t c_28 = C_f*(C_0[WCoef::C2]*conj(C_0[WCoef::C8]) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP8]));
	complex_t c_19 = C_tau_2*(C_0[WCoef::C1]*conj(C9_0) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP9]));
	complex_t c_29 = C_f*(C_0[WCoef::C2]*conj(C9_0) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP9]));

    complex_t w_22 = c_11 + c_12 + c_22;
    complex_t w_27 = c_17 + c_27;
    complex_t w_28 = c_18 + c_28;
    complex_t w_29 = c_19 + c_29;

    auto f = [&] (double w) -> double {
        complex_t D23 = Delta_i_23(s, cache.z, w);
        complex_t D27 = Delta_i_27(s, cache.z, w);
        return real(w_22 * tau_22(s, w, D23, D27) + 2.0 * (w_27 * tau_27(s, w, D23, D27) + w_28 * tau_28(s, w, D23, D27) + w_29 * tau_29(s, w, D23, D27)));
    };

    return integrate(f, s, 1, 1e-2);
}

double BXsllDecay::delta_bremB(double s) {
    if (s < 0.5) {
        return lerp(s, cache.delta_brems_lookup_low, cache.s_hat_low_bound.first, cache.s_hat_low_bound.second);
    } else {
        return lerp(s, cache.delta_brems_lookup_high, cache.s_hat_high_bound.first, cache.s_hat_high_bound.second);
    }
}

double BXsllDecay::delta_A_brem(double s) {
    auto C_0 = cache.C_LO;
    auto CP_0 = cache.C_LO;
    complex_t C10 = cache.C[WCoef::C10];
    complex_t CP10 = cache.C[WCoef::CP10];

    complex_t W_210 = (C_0[WCoef::C2] - C_0[WCoef::C1] / 6.) * C10 + (CP_0[WCoef::CP2] - CP_0[WCoef::CP1] / 6.) * CP10;
    complex_t W_810 = C_0[WCoef::C8] * C10 + CP_0[WCoef::CP8] * CP10;

    return pow(1 - s, 2) * real(W_810 * tau_810(s) + W_210 * tau_210(s, cache.z));
}

double BXsllDecay::delta_em(double s, double L_l) {
    auto C = cache.C;
    auto Cp = cache.C;

    double C_F = ObsQCDProxy().get_constants()->C_F;
    complex_t W_2 = pow(abs(C[WCoef::C2] + C_F * C[WCoef::C1]), 2) + pow(abs(Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]), 2);
    complex_t W_7 = pow(abs(C[WCoef::C7]), 2) + pow(abs(Cp[WCoef::CP7]), 2);
    complex_t W_9 = pow(abs(C[WCoef::C9]), 2) + pow(abs(Cp[WCoef::CP9]), 2);
    complex_t W_10 = pow(abs(C[WCoef::C10]), 2) + pow(abs(Cp[WCoef::CP10]), 2);
    complex_t W_27 = (C[WCoef::C2] + C_F * C[WCoef::C1]) * conj(C[WCoef::C7]) + (Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]) * conj(Cp[WCoef::CP7]);
    complex_t W_29 = (C[WCoef::C2] + C_F * C[WCoef::C1]) * conj(C[WCoef::C9]) + (Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]) * conj(Cp[WCoef::CP9]);
    complex_t W_79 = C[WCoef::C7] * conj(C[WCoef::C9]) + Cp[WCoef::CP7] * conj(Cp[WCoef::CP9]);

    return pow(1 - s, 2) * real(
        8 * (1 + 2 * s) * (
            W_9 * omega_99(s, L_l) +
            W_10 * omega_1010(s, L_l) + 
            real(W_29 * omega_29(s, L_l)) +
            W_2 * omega_22(s, L_l)
        ) +
        96 * real(
            W_79 * omega_79(s, L_l) +
            W_27 * omega_27(s, L_l)
        ) +
        8 * (4 + 8 / s) * W_7 * omega_77(s, L_l)
    );
}

double BXsllDecay::delta_A_em(double s, double L_l) {
    auto C = cache.C;
    auto Cp = cache.C;

    double C_F = ObsQCDProxy().get_constants()->C_F;
    complex_t W_210 = (C[WCoef::C2] + C_F * C[WCoef::C1]) * conj(C[WCoef::C10]) + (Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]) * conj(Cp[WCoef::CP10]);
    complex_t W_710 = C[WCoef::C7] * conj(C[WCoef::C10]) + Cp[WCoef::CP7] * conj(Cp[WCoef::CP10]);
    complex_t W_910 = C[WCoef::C9] * conj(C[WCoef::C10]) + Cp[WCoef::CP9] * conj(Cp[WCoef::CP10]);

    return pow(1 - s, 2) * real(
        -48. * W_710 * omega_710(s, L_l) +
        -24. * s * (
            W_910 * omega_910(s, L_l) +
            W_210 * omega_210(s, L_l, cache.z)
        )
    );
}

double BXsllDecay::A_FB(double s, double ml_hat, double L_l) {
    return -3 * A_FB_0(s, ml_hat) * (cache.pref_A0_0 + cache.pref_A0_1 * s / pow(1 - s, 2)) +
            cache.pref_delta_mb2 * delta_A_mb2(s) + 
            -3 / 8 * cache.pref_delta_mc2 * delta_A_mc2(s) +
            8 / 3 * cache.pref_delta_brems * delta_A_brem(s) +
            cache.pref_delta_em * delta_A_em(s, L_l);
}

double BXsllDecay::dB_ds(double s, double ml_hat, double L_l) {   
    return cache.pref_dB0_ds * dB0_ds(s, ml_hat) +
            cache.pref_delta_mb2 * delta_mb2(s) + 
            cache.pref_delta_mb3 * delta_mb3(s) + 
            cache.pref_delta_mc2 * delta_mc2(s) +
            cache.pref_delta_brems * (delta_bremA(s) + delta_bremB(s)) +
            cache.pref_delta_em * delta_em(s, L_l);
}

double BXsllDecay::BR_B_Xsll(double s_min, double s_max, int gen) {
    double ml_hat = gen == 2 ? cache.m_mu_hat : cache.m_tau_hat;
    double L_l = gen == 2 ? cache.L_l_mu : cache.L_l_tau;

    auto f = [&] (double s) {
        return dB_ds(s, ml_hat, L_l);
    };

    double res = cache.pref_dB_ds * integrate(f, s_min, s_max, 1e-3);
    return res;
}

// double BXsllDecay::BR_B_Xsll(double s_min, double s_max, int gen) {
//     double s_hat_min = s_min / pow(cache.m_b_1S, 2);
//     double s_hat_max = s_max / pow(cache.m_b_1S, 2);
//     double ml_hat = gen == 2 ? cache.m_mu_hat : cache.m_tau_hat;
//     double L_l = gen == 2 ? cache.L_l_mu : cache.L_l_tau;

//     std::ofstream fs;
//     fs.open("dB_ds_mu_high.csv", std::ios_base::app);
//     // fs << "s_hat,dB0_dS,delta_mb2,delta_mb3,delta_mc2,delta_brems_A,delta_brems_B,delta_em,delta_tot,W_7,W_9,W_10,W_79\n";

//     auto f = [&] (double s) {
//         double res = dB_ds(s, ml_hat, L_l);
//         fs << s << "," 
//            << cache.pref_dB_ds * cache.pref_dB0_ds * dB0_ds(s, ml_hat) << ","
//            << cache.pref_dB_ds * cache.pref_delta_mb2 * delta_mb2(s) << ","
//            << cache.pref_dB_ds * cache.pref_delta_mb3 * delta_mb3(s) << ","
//            << cache.pref_dB_ds * cache.pref_delta_mc2 * delta_mc2(s) << ","
//            << cache.pref_dB_ds * cache.pref_delta_brems * delta_bremA(s) << ","
//            << cache.pref_dB_ds * cache.pref_delta_brems * delta_bremB(s) << ","
//            << cache.pref_dB_ds * cache.pref_delta_em * delta_em(s, L_l) << ","
//            << cache.pref_dB_ds * (cache.pref_delta_mb2 * delta_mb2(s) + cache.pref_delta_mb3 * delta_mb3(s) + cache.pref_delta_mc2 * delta_mc2(s) + cache.pref_delta_brems * (delta_bremA(s) + delta_bremB(s)) + cache.pref_delta_em * delta_em(s, L_l)) << ","
//            << W_7(s) << ","
//            << W_9(s) << ","
//            << W_10(s) << ","
//            << W_79(s)
//            << "\n";
//         return res;
//     };

//     double res = cache.pref_dB_ds * integrate(f, s_hat_min, s_hat_max, 1e-3);
//     C7_new_base(0.19961, false);
    
//     return res;
// }


std::vector<ObservableValue> BXsllDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_B__Xs_mu_mu__LOW_Q2:   
        value = BR_B_Xsll(cache.s_hat_low_bound.first, cache.s_hat_low_bound.second, 2);
        break;
    case Observables::BR_B__Xs_mu_mu__HIGH_Q2:   
        value = BR_B_Xsll(cache.s_hat_high_bound.first, cache.s_hat_high_bound.second, 2);
        break;
    case Observables::BR_B__Xs_tau_tau__HIGH_Q2:   
        value = BR_B_Xsll(cache.s_hat_high_bound.first, cache.s_hat_high_bound.second, 3);
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BXsllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}