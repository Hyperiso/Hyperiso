#include "BXsllDecay.h"

double BXsllDecay::alpha_s(double mu) {
    return ObsQCDProxy()(AlphasConfig{mu, MassType::POLE, MassType::POLE});
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

double BXsllDecay::kappa(double f, double h, double alpha_s_mu_b) {
    return 1 - 2. * alpha_s_mu_b * h / (3. * PI * f);
}

double BXsllDecay::m_hat(double m, double mb) {
    return m / mb;
}

double BXsllDecay::z(double mc_hat) {
    return pow(mc_hat, 2);
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

scalar_t BXsllDecay::Gm1(double t) {
    scalar_t L = log((sqrt((complex_t)t)+sqrt((complex_t)t-4.))/2.);
    return -2.*I*PI*L-PI2/2.+2.*pow(L,2.);
}

scalar_t BXsllDecay::G0(double t) {
    scalar_t rt = sqrt(((complex_t)t-4.)/(complex_t)t);
    return -I*PI*rt-2.+2.*rt*log((sqrt((complex_t)t)+sqrt((complex_t)t-4.))/2.);
}

scalar_t BXsllDecay::Delta_i_23(double s, double z, double w) {
    return -2.+4./(w-s)*(z*Gm1(s/z)-z*Gm1(w/z)-s/2.*G0(s/z)+s/2.*G0(w/z));
}

scalar_t BXsllDecay::Delta_i_27(double s, double z, double w) {
    return 2.*(G0(s/z) - G0(w/z));
}

double BXsllDecay::tau_22(double s, double w, scalar_t Delta_23, scalar_t Delta_27) {
    return 8./27.*(w-s)*(1.-w)*(1.-w)/s/w/w/w*((3.*w*w+2.*s*s*(2.+w)-s*w*(5.-2.*w))*pow(abs(Delta_23),2.)
	+(2.*s*s*(2.+w)+s*w*(1.+2.*w))*pow(abs(Delta_27),2.)
	+4.*s*(w*(1.-w)-s*(2.+w))*real(Delta_23*conj(Delta_27)));
}

scalar_t BXsllDecay::tau_27(double s, double w, scalar_t Delta_23, scalar_t Delta_27) {
    return 8./3./s/w*(((1.-w)*(4.*s*s-s*w+w*w)+s*w*(4.+s-w)*log(w))*Delta_23
	-4.*s*s*(1.-w)+s*w*(4.+s-w)*log(w)*Delta_27);
}

scalar_t BXsllDecay::tau_28(double s, double w, scalar_t Delta_23, scalar_t Delta_27) {
    return 8./9./s/w/(w-s)*((pow(w-s,2.)*(2.*s-w)*(1.-w))*Delta_23
	-(2.*s*pow(w-s,2.)*(1.-w))*Delta_27
	+s*w*((1.+2.*s-2.*w)*Delta_23-2.*(1.+s-w)*Delta_27)*log(s/((1.+s-w)*(w*w+s*(1.-w)))));
}

scalar_t BXsllDecay::tau_29(double s, double w, scalar_t Delta_23, scalar_t Delta_27) {
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

scalar_t BXsllDecay::tau_210(double s, double z) {
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

scalar_t BXsllDecay::F(double r) {
    if(r < 1) {
        return 3./2./r*(1./sqrt(r*(1.-r))*atan(sqrt(r/(1.-r)))-1.);
    } 
    return 3./2./r*(1./2./sqrt(r*(r-1.))*(log((1.-sqrt(1.-1./r))/(1.+sqrt(1.-1./r)))+I*PI)-1.);
}

scalar_t BXsllDecay::Sigma_1(double s) {
    if (abs(s) < 0.4)
        return {
            23.787-120.948*s+365.373*s*s-584.206*s*s*s,
            1.653+6.009*s-17.080*s*s+115.880*s*s*s
        };

    scalar_t d = 1 - s;
    return {
        -148.061*d*d+492.539*d*d*d-1163.847*pow(d,4.)+1189.528*pow(d,5.)
        -261.287*d*d+1170.856*d*d*d-2546.948*pow(d,4.)+2540.023*pow(d,5.)
    };
}

double BXsllDecay::Sigma_2(double s) {
    if (abs(s) < 0.4)
        return 11.488-36.987*s+255.330*s*s-812.388*s*s*s+1011.791*s*s*s*s;

    scalar_t d = 1 - s;
    return -221.904*d*d+900.822*d*d*d-2031.620*pow(d,4.)+1984.303*pow(d,5.);
}

scalar_t BXsllDecay::Sigma_3(double s) {
    if (abs(s) < 0.4)
        return {
            109.311-846.039*s+2890.115*s*s-4179.072*s*s*s,
            4.606+17.650*s-53.244*s*s+348.069*s*s*s
        };

    scalar_t d = 1 - s;
    return {
        -298.730*d*d+828.0675*d*d*d-2217.6355*pow(d,4.)+2241.792*pow(d,5.),
        -528.759*d*d+2095.723*d*d*d-4681.843*pow(d,4.)+5036.677*pow(d,5.)
    };
}

scalar_t BXsllDecay::Sigma_7(double s, double z) {
    if (abs(s) < 0.4) {
        scalar_t a = pow(4 * z, 2); 
        return {
            -0.259023-28.424*s+205.533*s*s-603.219*s*s*s+722.031*s*s*s*s,
            (-12.20658-215.8208*(s-a)+412.1207*(s-a)*(s-a))*(s-a)*(s-a)*(s>a)
        };
    }
        
    scalar_t d = 1 - s;
    return {
        77.0256*d*d-264.705*d*d*d+595.814*pow(d,4.)-610.1637*pow(d,5.),
        135.858*d*d-618.990*d*d*d+1325.040*pow(d,4.)-1277.170*pow(d,5.)
    };
}

double BXsllDecay::omega_22(double s, double L_l, double L_b_5_GeV) {
    return L_l * (Sigma_2(s)/8./(1.-s)/(1.-s)/(1.+2.*s)+real(Sigma_1(s))/9./(1.-s)/(1.-s)/(1.+2.*s)*L_b_5_GeV)
	+64./81.*omega_1010(s, L_l) * L_b_5_GeV * L_b_5_GeV;
}

scalar_t BXsllDecay::omega_27(double s, double L_l, double L_b_5_GeV) {
    return L_l*(Sigma_3(s)/96./(1.-s)/(1.-s))+8./9.*omega_79(s, L_l)*L_b_5_GeV;
}

scalar_t BXsllDecay::omega_29(double s, double L_l, double L_b_5_GeV) {
    return L_l * (Sigma_1(s)/8./(1.-s)/(1.-s)/(1.+2.*s))+16./9.*omega_1010(s,L_l)*L_b_5_GeV;
}

scalar_t BXsllDecay::omega_210(double s, double L_l, double L_b_5_GeV, double z) {
    return L_l*(-Sigma_7(s, z)/24./s/(1.-s)/(1.-s))+8./9.*omega_910(s,L_l)*L_b_5_GeV;
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

scalar_t BXsllDecay::g(double z, double s) {
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

double BXsllDecay::breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had, double m_b) {
    double m_V_hat = m_V / m_b;
    double gamma_tot_hat = gamma_tot / m_b;
    double gamma_had_hat = gamma_had / m_b;
    return br * gamma_tot_hat * gamma_had_hat / (pow(s - pow(m_V_hat, 2.), 2.) + pow(m_V_hat * gamma_tot_hat, 2.));
}

double BXsllDecay::R_cc_cont(double s) {
    return s > 0.6 ? (s > 0.69 ? 1.02 : 11.33 * s - 6.8) : 0;
}

double BXsllDecay::R_cc(double s, double inv_alpha_em, double m_b) {
    double R_cc_res;

    for (size_t k = 0; k < cc_res_mass.size(); k++) {
        R_cc_res += breit_wigner(s, cc_res_mass[k], cc_res_br[k], cc_res_width_tot[k], cc_res_width_had[k], m_b);
    }

    return 9. * s * pow(inv_alpha_em, 2) * R_cc_res + R_cc_cont(s);
}

double BXsllDecay::PV_breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had, double m_b, double m_D_hat, double inv_alpha_em) {
    double s_c = 4 * pow(m_D_hat, 2);
    double m_V_hat = m_V / m_b;
    double m_V_hat2 = pow(m_V_hat, 2);
    double gamma_tot_hat = gamma_tot / m_b;
    double gamma_had_hat = gamma_had / m_b; 
    double den = (pow(s - pow(m_V_hat, 2.), 2.) + pow(m_V_hat * gamma_tot_hat, 2.));
    double B = breit_wigner(s, m_V, br, gamma_tot, gamma_had, m_b);
    return 9 * s * pow(inv_alpha_em, 2) * B * (0.5 * log(den / pow(s_c - s, 2)) + (s - m_V_hat2) / gamma_tot_hat * m_V_hat * (atan((s_c - m_V_hat2) / gamma_tot_hat * m_V_hat) - PI / 2));
}

double BXsllDecay::PV_R_cc_cont(double s, double m_D_hat) {
    double s_c = 4 * pow(m_D_hat, 2);
    return 1 / 3. * ((11.33 * s - 6.8) * log(abs((0.69 - s) / (s_c - s))) - 1.02 * log(abs(0.69 - s)) - 6.8 * log(s_c) - 2.902017); 
}

double BXsllDecay::PV_R_cc(double s, double inv_alpha_em, double m_b, double m_D_hat) {
    double PV_res;

    for (size_t k = 0; k < cc_res_mass.size(); k++) {
        PV_res += PV_breit_wigner(s, cc_res_mass[k], cc_res_br[k], cc_res_width_tot[k], cc_res_width_had[k], m_b, m_D_hat, inv_alpha_em);
    }

    return PV_res + PV_R_cc_cont(s, m_D_hat);
}

scalar_t BXsllDecay::g_ld(double z, double s, double inv_alpha_em, double m_D_hat, double m_b) {
    return g(z, 0) + (s * PV_R_cc(s, inv_alpha_em, m_b, m_D_hat) + I * PI * R_cc(s, inv_alpha_em, m_b)) / 3.;
}

scalar_t BXsllDecay::C9_eff_base(double s, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b, QCDOrder order, bool prime) {
    s = std::clamp(s, 1e-6, 1. - 1e-6);
    auto C = order == QCDOrder::LO ? (prime ? BPrime_FR_wilson_cache_LO : B_FR_wilson_cache_LO) : (prime ? BPrime_FR_wilson_cache : B_FR_wilson_cache);
    auto C_ids = WCoefMapper::get_group(prime ? WGroup::BPrime : WGroup::B);
    scalar_t g_0 = g(0, s);
    scalar_t g_1 = g(1, s);
    scalar_t g_mc = g_ld(mc_hat, s, inv_alpha_em, m_D_hat, m_b);

    return C[C_ids[8]]
	        +(-32./27.*C[C_ids[0]]-8./9.*C[C_ids[1]]-16./9.*C[C_ids[2]]+32./27.*C[C_ids[3]]-112./9.*C[C_ids[4]]+512./27.*C[C_ids[5]]) * L_mu
	        +4./3.*C[C_ids[2]]+64./9.*C[C_ids[4]]+64./27.*C[C_ids[5]]
	        +g_mc*(4./3.*C[C_ids[0]]+C[C_ids[1]]+6.*C[C_ids[2]]+60.*C[C_ids[4]])
	        +g_1*(-7./2.*C[C_ids[2]]-2./3.*C[C_ids[3]]-38.*C[C_ids[4]]-32./3.*C[C_ids[5]])
	        +g_0*(-1./2.*C[C_ids[2]]-2./3.*C[C_ids[3]]-8.*C[C_ids[4]]-32./3.*C[C_ids[5]]);
}

scalar_t BXsllDecay::C9_eff_LO(double s, bool prime) {
    return lerp(s, prime ? CP9_eff_LO_lookup : C9_eff_LO_lookup);
}

scalar_t BXsllDecay::F_17(double s) {
    return lerp(s, F_17_lookup);
}

scalar_t BXsllDecay::F_27(double s) {
    return lerp(s, F_27_lookup);
}

scalar_t BXsllDecay::F_19(double s) {
    return lerp(s, F_19_lookup);
}

scalar_t BXsllDecay::F_29(double s) {
    return lerp(s, F_29_lookup);
}

scalar_t BXsllDecay::C7_new_base(double s, double alpha_s_mu_b, double L_mu, bool prime) {
    auto C_0 = prime ? BPrime_FR_wilson_cache_LO : B_FR_wilson_cache_LO;
    auto C7_eff = prime ? BPrime_FR_wilson_cache[WCoef::CP7] : B_FR_wilson_cache[WCoef::C7];
    LOG_INFO("s_hat =", s);
    LOG_INFO("C7eff =", C7_eff);
    LOG_INFO("sigma_7 =", sigma_7(s, L_mu));
    LOG_INFO("F_17 =", F_17(s));
    LOG_INFO("F_27 =", F_27(s));
    LOG_INFO("F_87 =", f_87(s, L_mu));
    return (1.+alpha_s_mu_b/PI*sigma_7(s, L_mu))*C7_eff
	-alpha_s_mu_b/4./PI*(C_0.at(prime ? WCoef::CP1 :WCoef::C1)*F_17(s)+C_0.at(prime ? WCoef::CP2 :WCoef::C2)*F_27(s)+C_0.at(prime ? WCoef::CP8 :WCoef::C8)*f_87(s, L_mu));
}

scalar_t BXsllDecay::C9_new_base(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b, bool prime)
{
    if (abs(s - 1) < 1e-6) s = 1;
    auto C_0 = prime ? BPrime_FR_wilson_cache_LO : B_FR_wilson_cache_LO;
    return (1.+alpha_s_mu_b/PI*sigma_9(s))*C9_eff_base(s, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, this->w_config.order, prime)
	        -alpha_s_mu_b/4./PI*(C_0.at(prime ? WCoef::CP1 :WCoef::C1)*F_19(s)+C_0.at(prime ? WCoef::CP2 :WCoef::C2)*F_29(s)+C_0.at(prime ? WCoef::CP8 :WCoef::C8)*f_89(s));
}

scalar_t BXsllDecay::C10_new_base(double s, double alpha_s_mu_b, bool prime) {
    s = std::clamp(s, 1e-6, 1. - 1e-6);
    LOG_INFO("C10eff =", prime ? BPrime_FR_wilson_cache[WCoef::CP10] : B_FR_wilson_cache[WCoef::C10]);
    return (1.+alpha_s_mu_b/PI*sigma_9(s))*(prime ? BPrime_FR_wilson_cache[WCoef::CP10] : B_FR_wilson_cache[WCoef::C10]);
}

scalar_t BXsllDecay::C7_new(double s, bool prime) {
    return lerp(s, prime ? CP7_new_lookup : C7_new_lookup);
}

scalar_t BXsllDecay::C9_new(double s, bool prime) {
    return lerp(s, prime ? CP9_new_lookup : C9_new_lookup);
}

scalar_t BXsllDecay::C10_new(double s, bool prime) {
    return lerp(s, prime ? CP10_new_lookup : C10_new_lookup);
}

double BXsllDecay::W_7(double s, double alpha_s_mu_b, double L_mu) {
    return pow(abs(C7_new_base(s, alpha_s_mu_b, L_mu, false)), 2) + pow(abs(C7_new_base(s, alpha_s_mu_b, L_mu, true)), 2);
}

double BXsllDecay::W_9(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    return pow(abs(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, false)), 2) + pow(abs(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, true)), 2);
}

double BXsllDecay::W_10(double s, double alpha_s_mu_b) {
    return pow(abs(C10_new_base(s, alpha_s_mu_b, false)), 2) + pow(abs(C10_new_base(s, alpha_s_mu_b, true)), 2);
}

scalar_t BXsllDecay::W_27(double s, double alpha_s_mu_b, double L_mu) {
    return B_FR_wilson_cache[WCoef::C2] * conj(C7_new_base(s, alpha_s_mu_b, L_mu, false)) 
            + BPrime_FR_wilson_cache[WCoef::CP2] * conj(C7_new_base(s, alpha_s_mu_b, L_mu, true));
}

scalar_t BXsllDecay::W_29(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    return B_FR_wilson_cache[WCoef::C2] * conj(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, false)) +
            BPrime_FR_wilson_cache[WCoef::CP2] * conj(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, true));
}

scalar_t BXsllDecay::W_210(double s, double alpha_s_mu_b) {
    return B_FR_wilson_cache[WCoef::C2] * conj(C10_new_base(s, alpha_s_mu_b, false)) +
            BPrime_FR_wilson_cache[WCoef::CP2] * conj(C10_new_base(s, alpha_s_mu_b, true));
}

double BXsllDecay::W_79(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    return real(C7_new_base(s, alpha_s_mu_b, L_mu, false) * conj(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, false))) 
            + real(C7_new_base(s, alpha_s_mu_b, L_mu, true) * conj(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, true)));
}

double BXsllDecay::W_710(double s, double alpha_s_mu_b, double L_mu) {
    return real(C7_new_base(s, alpha_s_mu_b, L_mu, false) * conj(C10_new_base(s, alpha_s_mu_b, false))) 
            + real(C7_new_base(s, alpha_s_mu_b, L_mu, true) * conj(C10_new_base(s, alpha_s_mu_b, true)));
}

double BXsllDecay::W_910(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    return real(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, false) * conj(C10_new_base(s, alpha_s_mu_b, false))) 
            + real(C9_new_base(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b, true) * conj(C10_new_base(s, alpha_s_mu_b, true)));
}

scalar_t BXsllDecay::H7(double s, double ml_hat, double alpha_s_mu_b) {
    return 4 * (1. + 2. * ml_hat * ml_hat / s) * (1. + 2. / s) * (1. + alpha_s_mu_b / PI * tau_77(s));
}

scalar_t BXsllDecay::H9(double s, double ml_hat, double alpha_s_mu_b) {
    return (1. + 2. * ml_hat * ml_hat / s) * (1. + 2. * s) * (1. + alpha_s_mu_b / PI * tau_99(s));
}

scalar_t BXsllDecay::H10(double s, double ml_hat, double alpha_s_mu_b) {
    return ((1.+2.*s)+2.*ml_hat*ml_hat/s*(1.-4.*s))*(1.+alpha_s_mu_b/PI*tau_99(s));
}

scalar_t BXsllDecay::H79(double s, double ml_hat, double alpha_s_mu_b) {
    return 12 * (1.+2.*ml_hat*ml_hat/s)*(1.+alpha_s_mu_b/PI*tau_79(s));
}

double BXsllDecay::dB0_ds(double s, double ml_hat, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    double W_Q1 = pow(abs(BScalar_FR_wilson_cache[WCoef::CQ1]), 2) + pow(abs(BPrime_FR_wilson_cache[WCoef::CPQ1]), 2);
    double W_Q2 = pow(abs(BScalar_FR_wilson_cache[WCoef::CQ2]), 2) + pow(abs(BPrime_FR_wilson_cache[WCoef::CPQ2]), 2);
    double W_10Q2 = real(BScalar_FR_wilson_cache[WCoef::CQ2] * conj(C10_new_base(s, alpha_s_mu_b, false))) + real(BPrime_FR_wilson_cache[WCoef::CPQ2] * conj(C10_new_base(s, alpha_s_mu_b, true)));

    return pow(1 - s, 2) * sqrt(1 - 4 * pow(ml_hat, 2) / s) * (
            H7(s, ml_hat, alpha_s_mu_b) * W_7(s, alpha_s_mu_b, L_mu) + 
            H9(s, ml_hat, alpha_s_mu_b) * W_9(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + 
            H10(s, ml_hat, alpha_s_mu_b) * W_10(s, alpha_s_mu_b) + 
            H79(s, ml_hat, alpha_s_mu_b) * W_79(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) +
            1.5 * (s - 4 * ml_hat * ml_hat) * W_Q1 +
            1.5 * s * W_Q2 +
            6 * ml_hat * W_10Q2
        );
}

double BXsllDecay::A_FB_0(double s, double ml_hat, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    double W_7Q1 = real(C7_new(s, false) * BScalar_FR_wilson_cache[WCoef::CQ1]) + real(C7_new(s, true) * BPrime_FR_wilson_cache[WCoef::CPQ1]);
    double W_9Q1 = real(C9_new(s, false) * BScalar_FR_wilson_cache[WCoef::CQ1]) + real(C9_new(s, true) * BPrime_FR_wilson_cache[WCoef::CPQ1]);

    return pow(1 - s, 2) * sqrt(1 - 4 * pow(ml_hat, 2) / s) * (
            2 * (1 + alpha_s_mu_b * tau_710(s) / PI) * W_710(s, alpha_s_mu_b, L_mu) + 
            s * (1 + alpha_s_mu_b * tau_910(s) / PI) * W_910(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + 
            ml_hat * (2 * W_7Q1 + W_9Q1)
        );
}

double BXsllDecay::delta_A_mb2(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    return s * (9 + 14 * s - 15 * pow(s, 2)) * W_910(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + 2 * (7 + 10 * s - 9 * s * s) * W_710(s, alpha_s_mu_b, L_mu);
}

double BXsllDecay::delta_mb2(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    return -4 * (6 + 3 * s - 5 * pow(s, 3)) * W_7(s, alpha_s_mu_b, L_mu) / s + (1 - 15 * s * s + 10 * pow(s, 3)) * (W_9(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + W_10(s, alpha_s_mu_b)) - 4 * (5 + 6 * s - 7 * s * s) * W_79(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b);
}

double BXsllDecay::delta_mb3(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    if (s > 0.4)
        return 0;

    return (5.*pow(s,4.)+19.*pow(s,3.)+9.*s*s-7.*s+22.)/6./(1.-s)*4.*W_7(s, alpha_s_mu_b, L_mu)/s+
            (10.*pow(s,4.)+23.*pow(s,3.)-9.*s*s+13.*s+11.)/6./(1.-s)*(W_9(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + W_10(s, alpha_s_mu_b))+
            4.*(-3.*pow(s,3.)+17.*s*s-s+3.)/2./(1.-s)* W_79(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b);
}

double BXsllDecay::delta_mc2(double s, double z, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b) {
    scalar_t f = F(s / (4. * z));
    return pow(1 - s, 2) * real((1 + 6 * s - s * s) * f * W_27(s, alpha_s_mu_b, L_mu) / s + (2 + s) * f * W_29(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b));
}

double BXsllDecay::delta_A_mc2(double s, double z, double alpha_s_mu_b) {
    return pow(1 - s, 2) * real((1 + 3 * s) * F(s / (4. * z)) * W_210(s, alpha_s_mu_b));
}

double BXsllDecay::delta_bremA(double s) {
    scalar_t C7_0 = B_FR_wilson_cache_LO[WCoef::C7];
    scalar_t C8_0 = B_FR_wilson_cache_LO[WCoef::C8];
    scalar_t CP7_0 = BPrime_FR_wilson_cache_LO[WCoef::CP7];
    scalar_t CP8_0 = BPrime_FR_wilson_cache_LO[WCoef::CP8];
    scalar_t C9_0 = C9_eff_LO(s, false);
    scalar_t CP9_0 = C9_eff_LO(s, true);

    scalar_t c_78 = ObsQCDProxy().get_constants()->C_F * (C7_0 * conj(C8_0) + CP7_0 * conj(CP8_0));
    scalar_t c_88 = ObsQCDProxy().get_constants()->C_F * (C8_0 * conj(C8_0) + CP8_0 * conj(CP8_0));
    scalar_t c_89 = ObsQCDProxy().get_constants()->C_F * (C8_0 * conj(C9_0) + CP8_0 * conj(CP9_0));

    return 2 * real(c_78 * tau_78(s) + c_89 * tau_89(s)) + c_88 * tau_88(s);
}

double BXsllDecay::delta_bremB_base(double t, double z, double s_min, double s_max) {
    double s = t * s_max + (1. - t) * s_min;
    auto C_0 = B_FR_wilson_cache_LO;
    auto CP_0 = BPrime_FR_wilson_cache_LO;
    scalar_t C9_0 = C9_eff_LO(s, false);
    scalar_t CP9_0 = C9_eff_LO(s, true);

    double C_f = ObsQCDProxy().get_constants()->C_F;
    double C_tau_1 = C_f / (4 * pow(ObsQCDProxy().get_constants()->Nc, 2));
    double C_tau_2 = -C_f / (2 * ObsQCDProxy().get_constants()->Nc);

    scalar_t c_11 = C_tau_1 * (C_0[WCoef::C1] * conj(C_0[WCoef::C1]) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP1]));
	scalar_t c_12 = 2 * C_tau_2 * real((C_0[WCoef::C1]*conj(C_0[WCoef::C2])) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP2]));
	scalar_t c_22 = C_f*(C_0[WCoef::C2]*conj(C_0[WCoef::C2]) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP2]));
	scalar_t c_17 = C_tau_2*(C_0[WCoef::C1]*conj(C_0[WCoef::C7]) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP7]));
	scalar_t c_27 = C_f*(C_0[WCoef::C2]*conj(C_0[WCoef::C7]) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP7]));
	scalar_t c_18 = C_tau_2*(C_0[WCoef::C1]*conj(C_0[WCoef::C8]) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP8]));
	scalar_t c_28 = C_f*(C_0[WCoef::C2]*conj(C_0[WCoef::C8]) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP8]));
	scalar_t c_19 = C_tau_2*(C_0[WCoef::C1]*conj(C9_0) + CP_0[WCoef::CP1] * conj(CP_0[WCoef::CP9]));
	scalar_t c_29 = C_f*(C_0[WCoef::C2]*conj(C9_0) + CP_0[WCoef::CP2] * conj(CP_0[WCoef::CP9]));

    scalar_t w_22 = c_11 + c_12 + c_22;
    scalar_t w_27 = c_17 + c_27;
    scalar_t w_28 = c_18 + c_28;
    scalar_t w_29 = c_19 + c_29;

    auto f = [&] (double w) {
        scalar_t D23 = Delta_i_23(s, z, w);
        scalar_t D27 = Delta_i_27(s, z, w);
        return w_22 * tau_22(s, w, D23, D27) + 2 * real(w_27 * tau_27(s, w, D23, D27) + w_28 * tau_28(s, w, D23, D27) + w_29 * tau_29(s, w, D23, D27));
    };

    return integrate(f, s, 1, 1e-2);
}

double BXsllDecay::delta_bremB(double s, double m_b) {
    if (s < 0.5) {
        double t = (s * pow(m_b, 2) - q2_low_bound.second) / (q2_low_bound.first - q2_low_bound.second);
        return lerp(t, delta_brems_lookup_low);
    } else {
        double t = (s * pow(m_b, 2) - q2_high_bound.second) / (q2_high_bound.first - q2_high_bound.second);
        return lerp(t, delta_brems_lookup_high);
    }
}

double BXsllDecay::delta_A_brem(double s, double z) {
    auto C_0 = B_FR_wilson_cache_LO;
    auto CP_0 = BPrime_FR_wilson_cache_LO;
    scalar_t C10 = B_FR_wilson_cache[WCoef::C10];
    scalar_t CP10 = BPrime_FR_wilson_cache[WCoef::CP10];

    scalar_t W_210 = (C_0[WCoef::C2] - C_0[WCoef::C1] / 6.) * C10 + (CP_0[WCoef::CP2] - CP_0[WCoef::CP1] / 6.) * CP10;
    scalar_t W_810 = C_0[WCoef::C8] * C10 + CP_0[WCoef::CP8] * CP10;

    return pow(1 - s, 2) * real(W_810 * tau_810(s) + W_210 * tau_210(s, z));
}

double BXsllDecay::delta_em(double s, double L_l, double L_b_5_GeV) {
    auto C = B_FR_wilson_cache;
    auto Cp = BPrime_FR_wilson_cache;

    double C_F = ObsQCDProxy().get_constants()->C_F;
    scalar_t W_2 = pow(abs(C[WCoef::C2] + C_F * C[WCoef::C1]), 2) + pow(abs(Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]), 2);
    scalar_t W_7 = pow(abs(C[WCoef::C7]), 2) + pow(abs(Cp[WCoef::CP7]), 2);
    scalar_t W_9 = pow(abs(C[WCoef::C9]), 2) + pow(abs(Cp[WCoef::CP9]), 2);
    scalar_t W_10 = pow(abs(C[WCoef::C10]), 2) + pow(abs(Cp[WCoef::CP10]), 2);
    scalar_t W_27 = (C[WCoef::C2] + C_F * C[WCoef::C1]) * conj(C[WCoef::C7]) + (Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]) * conj(Cp[WCoef::CP7]);
    scalar_t W_29 = (C[WCoef::C2] + C_F * C[WCoef::C1]) * conj(C[WCoef::C9]) + (Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]) * conj(Cp[WCoef::CP9]);
    scalar_t W_79 = C[WCoef::C7] * conj(C[WCoef::C9]) + Cp[WCoef::CP7] * conj(Cp[WCoef::CP9]);

    return pow(1 - s, 2) * (
        8 * (1 + 2 * s) * (
            W_9 * omega_99(s, L_l) +
            W_10 * omega_1010(s, L_l) + 
            real(W_29 * omega_29(s, L_l, L_b_5_GeV)) +
            W_2 * omega_22(s, L_l, L_b_5_GeV)
        ) +
        96 * real(
            W_79 * omega_79(s, L_l) +
            W_27 * omega_27(s, L_l, L_b_5_GeV)
        ) +
        8 * (4 + 8 / s) * W_7 * omega_77(s, L_l)
    );
}

double BXsllDecay::delta_A_em(double s, double L_l, double L_b_5_GeV, double z) {
    auto C = B_FR_wilson_cache;
    auto Cp = BPrime_FR_wilson_cache;

    double C_F = ObsQCDProxy().get_constants()->C_F;
    scalar_t W_210 = (C[WCoef::C2] + C_F * C[WCoef::C1]) * conj(C[WCoef::C10]) + (Cp[WCoef::CP2] + C_F * Cp[WCoef::CP1]) * conj(Cp[WCoef::CP10]);
    scalar_t W_710 = C[WCoef::C7] * conj(C[WCoef::C10]) + Cp[WCoef::CP7] * conj(Cp[WCoef::CP10]);
    scalar_t W_910 = C[WCoef::C9] * conj(C[WCoef::C10]) + Cp[WCoef::CP9] * conj(Cp[WCoef::CP10]);

    return pow(1 - s, 2) * (
        -48 * W_710 * omega_710(s, L_l) +
        -24 * s * (
            W_910 * omega_910(s, L_l) +
            W_210 * omega_210(s, L_l, L_b_5_GeV, z)
        )
    );
}

double BXsllDecay::pref_A0_0(double lambda_2, double g_lam_z, double f_z, double m_b) {
    return 1 + 3 * lambda_2 * g_lam_z / (2 * pow(m_b, 2) * f_z);
}

double BXsllDecay::pref_A0_1(double lambda_1, double m_b) {
    return 4 * lambda_1 / (3 * pow(m_b, 2));
}

double BXsllDecay::A_FB(double s, double ml_hat, double alpha_s_mu_b, double z, double L_l, double L_b_5_GeV, double pref_A0_0, 
                double pref_A0_1, double pref_delta_mb2, double pref_delta_mc2, double pref_delta_brems, double pref_delta_em,
                double m_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat) 
{
    return -3 * A_FB_0(s, ml_hat, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) * (pref_A0_0 + pref_A0_1 * s / pow(1 - s, 2)) +
            pref_delta_mb2 * delta_A_mb2(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + 
            -3 / 8 * pref_delta_mc2 * delta_A_mc2(s, z, alpha_s_mu_b) +
            8 / 3 * pref_delta_brems * delta_A_brem(s, z) +
            pref_delta_em * delta_A_em(s, L_l, L_b_5_GeV, z);
}

double BXsllDecay::pref_dB0_ds(double lambda_2,
                                 double g_lam_z,
                                 double f_z,
                                 double m_b,
                                 double rho_1,
                                 double g_rho_z) {
    return 1. + 3. * lambda_2 * g_lam_z / (2. * pow(m_b, 2) * f_z) - rho_1 * g_rho_z / (6. * pow(m_b, 3) * f_z);
}

double BXsllDecay::pref_delta_mb2(double lambda_2,
                                    double m_b,
                                    scalar_t V_tb)
{
    return 3. * lambda_2 / (2. * pow(m_b * abs(V_tb), 2));
}

double BXsllDecay::pref_delta_mb3(double rho_1, double m_b, scalar_t V_tb) {
    return -rho_1 / (pow(m_b, 3) * pow(abs(V_tb), 2));
}

double BXsllDecay::pref_delta_mc2(double lambda_2,
                                    double m_c,
                                    scalar_t V_tb,
                                    scalar_t V_ts,
                                    scalar_t V_cb,
                                    scalar_t V_cs)
{
    return 8. * lambda_2 / (9. * pow(m_c, 2)) * abs(conj(V_cs) * V_cb / (conj(V_ts) * pow(V_tb, 3)));
}

double BXsllDecay::pref_delta_brems(double alpha_s_mu_b) {
    return alpha_s_mu_b / (4. * PI);
}

double BXsllDecay::pref_delta_em(double inv_alpha_em) {
    return 1. / (4. * PI * inv_alpha_em);
}


double BXsllDecay::dB_ds(double s, double ml_hat, double alpha_s_mu_b, double z,
                           double L_l, double L_b_5_GeV, double pref_dB0_ds, double pref_delta_mb2,
                           double pref_delta_mb3, double pref_delta_mc2, double pref_delta_brems, double pref_delta_em,
                           double m_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat)
{   
    return pref_dB0_ds * dB0_ds(s, ml_hat, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) +
            pref_delta_mb2 * delta_mb2(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + 
            pref_delta_mb3 * delta_mb3(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + 
            pref_delta_mc2 * delta_mc2(s, z, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) +
            pref_delta_brems * (delta_bremA(s) + delta_bremB(s, m_b)) +
            pref_delta_em * delta_em(s, L_l, L_b_5_GeV);
}

double BXsllDecay::ckm(scalar_t V_tb, scalar_t V_ts, scalar_t V_cb) {
    return pow(abs(V_tb * conj(V_ts) / V_cb), 2);
}

double BXsllDecay::pref(double BR_BXclnu, double inv_alpha_em, double f, double kappa, double ckm) {
    return BR_BXclnu * ckm / (pow(2 * PI * inv_alpha_em, 2) * f * kappa);
}

double BXsllDecay::BR_B_Xsll(double s_min, double s_max, double m_b, double pref, double ml_hat, double alpha_s_mu_b, 
                        double z, double L_l, double L_b_5_GeV, double pref_dB0_ds, double pref_delta_mb2,
                        double pref_delta_mb3, double pref_delta_mc2, double pref_delta_brems, double pref_delta_em, 
                        double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat)
{
    double s_hat_min = s_min / pow(m_b, 2);
    double s_hat_max = s_max / pow(m_b, 2);

    std::ofstream fs;
    fs.open("dB_ds_mu_high.csv", std::ios_base::app);
    // fs << "s_hat,dB0_dS,delta_mb2,delta_mb3,delta_mc2,delta_brems_A,delta_brems_B,delta_em,delta_tot,W_7,W_9,W_10,W_79\n";

    auto f = [&] (double s) {
        double res = dB_ds(s, ml_hat, alpha_s_mu_b, z, L_l, L_b_5_GeV, pref_dB0_ds, pref_delta_mb2, pref_delta_mb3, pref_delta_mc2, pref_delta_brems, pref_delta_em, m_b,  L_mu, mc_hat, inv_alpha_em, m_D_hat);
        fs << s << "," 
           << pref * pref_dB0_ds * dB0_ds(s, ml_hat, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) << ","
           << pref * pref_delta_mb2 * delta_mb2(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) << ","
           << pref * pref_delta_mb3 * delta_mb3(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) << ","
           << pref * pref_delta_mc2 * delta_mc2(s, z, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) << ","
           << pref * pref_delta_brems * delta_bremA(s) << ","
           << pref * pref_delta_brems * delta_bremB(s, m_b) << ","
           << pref * pref_delta_em * delta_em(s, L_l, L_b_5_GeV) << ","
           << pref * (pref_delta_mb2 * delta_mb2(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + pref_delta_mb3 * delta_mb3(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + pref_delta_mc2 * delta_mc2(s, z, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) + pref_delta_brems * (delta_bremA(s) + delta_bremB(s, m_b)) + pref_delta_em * delta_em(s, L_l, L_b_5_GeV)) << ","
           << W_7(s, alpha_s_mu_b, L_mu) << ","
           << W_9(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b) << ","
           << W_10(s, alpha_s_mu_b) << ","
           << W_79(s, alpha_s_mu_b, L_mu, mc_hat, inv_alpha_em, m_D_hat, m_b)
           << "\n";
        return res;
    };

    double res = pref * integrate(f, s_hat_min, s_hat_max, 1e-3);
    C7_new_base(0.19961, alpha_s_mu_b, L_mu, false);
    
    return res;
}

void BXsllDecay::build_op_tree() {
    // TODO: Make it possible to choose lepton gen

    // SM Parameters
    auto inv_alpha_em = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 1));
    auto m_mu = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 13));
    auto m_tau = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 15));
    auto m_c = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 4));
    auto m_b_1S = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "QCD", LhaID(5, 3)));

    auto V_tb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(2, 2)));
    auto V_ts = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(2, 1)));
    auto V_cb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(1, 2)));
    auto V_cs = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(1, 1)));

    // Flavor parameters
    auto m_D = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 411));

    // Wilson node
    auto wilson = this->get_wilson_node();
    auto mu_b = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, "B_SCALE", 1));

    // Misc experimental input
    auto BR_B_Xclnu = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Xs", 2));
    auto lambda_1 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Xsll", 1)); 
    auto lambda_2 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Xs", 6)); 
    auto rho_1 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Xsll", 2)); 
    
    // Operator nodes
    auto m_mu_hat = std::make_shared<OperatorNode>("m_mu_hat", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return m_hat(values[0], values[1]); });
    m_mu_hat->addChildren({m_mu, m_b_1S});
    auto m_tau_hat = std::make_shared<OperatorNode>("m_tau_hat", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return m_hat(values[0], values[1]); });
    m_tau_hat->addChildren({m_tau, m_b_1S});
    auto m_c_hat = std::make_shared<OperatorNode>("m_c_hat", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return m_hat(values[0], values[1]); });
    m_c_hat->addChildren({m_c, m_b_1S});
    auto m_D_hat = std::make_shared<OperatorNode>("m_D_hat", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return m_hat(values[0], values[1]); });
    m_D_hat->addChildren({m_D, m_b_1S});
    auto z = std::make_shared<OperatorNode>("z", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pow(values[0], 2); });
    z->addChildren({m_c_hat});
    auto L_b = std::make_shared<OperatorNode>("L_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return log(values[0] / values[1]); });
    L_b->addChildren({mu_b, m_b_1S});
    auto L_b_5GeV = std::make_shared<OperatorNode>("L_b_5GeV", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return log(values[0] / 5.); });
    L_b_5GeV->addChildren({mu_b});
    auto L_l_mu = std::make_shared<OperatorNode>("L_l_mu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return -2. * log(values[0]); });
    L_l_mu->addChildren({m_mu_hat});
    auto L_l_tau = std::make_shared<OperatorNode>("L_l_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return -2. * log(values[0]); });
    L_l_tau->addChildren({m_tau_hat});
    auto alpha_s_mu_b = std::make_shared<OperatorNode>("alpha_s(mu_b)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return alpha_s(values[0]); });
    alpha_s_mu_b->addChildren({mu_b});
    auto f_z = std::make_shared<OperatorNode>("f", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return f(values[0]); });
    f_z->addChildren({z});
    auto h_z = std::make_shared<OperatorNode>("h", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return h(values[0]); });
    h_z->addChildren({z});
    auto kappa_z = std::make_shared<OperatorNode>("kappa", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return kappa(values[0], values[1], values[2]); });
    kappa_z->addChildren({f_z, h_z, alpha_s_mu_b});
    auto g_lambda_z = std::make_shared<OperatorNode>("g_lambda", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return g_lambda(values[0]); });
    g_lambda_z->addChildren({z});
    auto g_rho_z = std::make_shared<OperatorNode>("g_rho", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return g_rho(values[0]); });
    g_rho_z->addChildren({z});

    auto wilson_cache = std::make_shared<OperatorNode>("wilson_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        BXsllDecay::B_FR_wilson_cache = this->w_proxy->getAFR(WGroup::B, this->w_config.order);
        BXsllDecay::BPrime_FR_wilson_cache = this->w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
        BXsllDecay::B_FR_wilson_cache_LO = this->w_proxy->getAFR(WGroup::B, QCDOrder::LO);
        BXsllDecay::BPrime_FR_wilson_cache_LO = this->w_proxy->getAFR(WGroup::BPrime, QCDOrder::LO);
        BXsllDecay::BScalar_FR_wilson_cache = this->w_proxy->getAFR(WGroup::BScalar, this->w_config.order);
        return 0; 
    });
    wilson_cache->addChildren({wilson});

    auto C9_eff_lookup = std::make_shared<OperatorNode>("C9_eff_LO_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BXsllDecay::C9_eff_base, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
        fill_cache(bound_func, 0, 1, C9_eff_LO_lookup, values[0], values[1], values[2], values[3], values[4], QCDOrder::LO, false); 
        fill_cache(bound_func, 0, 1, CP9_eff_LO_lookup, values[0], values[1], values[2], values[3], values[4], QCDOrder::LO, true); 
        return 0; 
    });
    C9_eff_lookup->addChildren({L_b, m_c_hat, inv_alpha_em, m_D_hat, m_b_1S, wilson_cache});

    auto n_f_17_lookup = std::make_shared<OperatorNode>("f_17_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        fill_cache(f_17, 0, 1, F_17_lookup, 0., 0.0625, FF_ORDER);
        return 0; 
    });
    n_f_17_lookup->addChildren({L_b, z});

    auto n_f_27_lookup = std::make_shared<OperatorNode>("f_27_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {
        fill_cache(f_27, 0, 1, F_27_lookup, 0., 0.0625, FF_ORDER);
        return 0; 
    });
    n_f_27_lookup->addChildren({L_b, z});

    auto n_f_19_lookup = std::make_shared<OperatorNode>("f_19_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        fill_cache(f_19_pole, 0, 1, F_19_lookup, 0., 0.0625, FF_ORDER);
        return 0; 
    });
    n_f_19_lookup->addChildren({L_b, z});

    auto n_f_29_lookup = std::make_shared<OperatorNode>("f_29_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        fill_cache(f_29_pole, 0, 1, F_29_lookup, 0., 0.0625, FF_ORDER);
        return 0; 
    });
    n_f_29_lookup->addChildren({L_b, z});

    // auto n_C7_new_lookup = std::make_shared<OperatorNode>("C7_new_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
    //     auto bound_func = std::bind(&BXsllDecay::C7_new_base, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    //     fill_cache(bound_func, C7_new_lookup, values[0], values[1], false); 
    //     fill_cache(bound_func, CP7_new_lookup, values[0], values[1], true); 
    //     return 0; 
    // });
    // n_C7_new_lookup->addChildren({L_b, m_c_hat, inv_alpha_em, m_D_hat, n_f_17_lookup, n_f_27_lookup, wilson_cache});

    // auto n_C9_new_lookup = std::make_shared<OperatorNode>("C9_new_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
    //     auto bound_func = std::bind(&BXsllDecay::C9_new_base, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
    //     fill_cache(bound_func, C9_new_lookup, values[0], values[1], values[2], values[3], values[4], values[5], false); 
    //     fill_cache(bound_func, CP9_new_lookup, values[0], values[1], values[2], values[3], values[4], values[5], true); 
    //     return 0; 
    // });
    // n_C9_new_lookup->addChildren({alpha_s_mu_b, L_b, m_c_hat, inv_alpha_em, m_D_hat, m_b_1S, n_f_19_lookup, n_f_29_lookup, wilson_cache});

    // auto n_C10_new_lookup = std::make_shared<OperatorNode>("C10_new_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
    //     auto bound_func = std::bind(&BXsllDecay::C10_new_base, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    //     fill_cache(bound_func, C10_new_lookup, values[0], false); 
    //     fill_cache(bound_func, CP10_new_lookup, values[0], true); 
    //     return 0; 
    // });
    // n_C10_new_lookup->addChildren({alpha_s_mu_b, wilson_cache});

    auto n_delta_brems_B_lookup = std::make_shared<OperatorNode>("delta_brems_B_lookup", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BXsllDecay::delta_bremB_base, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        fill_cache(bound_func, 0, 1, delta_brems_lookup_low, values[0], q2_low_bound.first / pow(values[1], 2), q2_low_bound.second / pow(values[1], 2));
        fill_cache(bound_func, 0, 1, delta_brems_lookup_high, values[0], q2_high_bound.first / pow(values[1], 2), q2_high_bound.second / pow(values[1], 2)); 
        return 0; 
    });
    n_delta_brems_B_lookup->addChildren({z, m_b_1S, wilson_cache});
    
    auto n_ckm = std::make_shared<OperatorNode>("ckm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return ckm(values[0], values[1], values[2]); });
    n_ckm->addChildren({V_tb, V_ts, V_cb});
    auto n_pref = std::make_shared<OperatorNode>("pref", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref(values[0], values[1], values[2], values[3], values[4]); });
    n_pref->addChildren({BR_B_Xclnu, inv_alpha_em, f_z, kappa_z, n_ckm});
    auto n_pref_dB0_ds = std::make_shared<OperatorNode>("pref_dB0_ds", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref_dB0_ds(values[0], values[1], values[2], values[3], values[4], values[5]); });
    n_pref_dB0_ds->addChildren({lambda_2, g_lambda_z, f_z, m_b_1S, rho_1, g_rho_z});
    auto n_pref_delta_mb2 = std::make_shared<OperatorNode>("pref_delta_mb2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref_delta_mb2(values[0], values[1], values[2]); });
    n_pref_delta_mb2->addChildren({lambda_2, m_b_1S, V_tb});
    auto n_pref_delta_mb3 = std::make_shared<OperatorNode>("pref_delta_mb3", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref_delta_mb3(values[0], values[1], values[2]); });
    n_pref_delta_mb3->addChildren({rho_1, m_b_1S, V_tb});
    auto n_pref_delta_mc2 = std::make_shared<OperatorNode>("pref_delta_mc2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref_delta_mc2(values[0], values[1], values[2], values[3], values[4], values[5]); });
    n_pref_delta_mc2->addChildren({lambda_2, m_c, V_tb, V_ts, V_cb, V_cs});
    auto n_pref_delta_brems = std::make_shared<OperatorNode>("pref_delta_brems", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref_delta_brems(values[0]); });
    n_pref_delta_brems->addChildren({alpha_s_mu_b});
    auto n_pref_delta_em = std::make_shared<OperatorNode>("pref_delta_em", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref_delta_em(values[0]); });
    n_pref_delta_em->addChildren({inv_alpha_em});
    auto br_low_mu = std::make_shared<OperatorNode>("BR(B>X_sll)_low_q2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return BR_B_Xsll(q2_low_bound.first, q2_low_bound.second, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16]); });
    br_low_mu->addChildren({m_b_1S, n_pref, m_mu_hat, alpha_s_mu_b, z, L_l_mu, L_b_5GeV, n_pref_dB0_ds, n_pref_delta_mb2, n_pref_delta_mb3, n_pref_delta_mc2, n_pref_delta_brems, n_pref_delta_em, L_b, m_c_hat, inv_alpha_em, m_D_hat, C9_eff_lookup, wilson_cache, n_delta_brems_B_lookup, n_f_17_lookup, n_f_19_lookup, n_f_27_lookup, n_f_29_lookup});
    auto br_high_mu = std::make_shared<OperatorNode>("BR(B>X_sll)_high_q2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return BR_B_Xsll(q2_high_bound.first, q2_high_bound.second, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16]); });
    br_high_mu->addChildren({m_b_1S, n_pref, m_mu_hat, alpha_s_mu_b, z, L_l_mu, L_b_5GeV, n_pref_dB0_ds, n_pref_delta_mb2, n_pref_delta_mb3, n_pref_delta_mc2, n_pref_delta_brems, n_pref_delta_em, L_b, m_c_hat, inv_alpha_em, m_D_hat, C9_eff_lookup, wilson_cache, n_delta_brems_B_lookup, n_f_17_lookup, n_f_19_lookup, n_f_27_lookup, n_f_29_lookup});
    auto br_high_tau = std::make_shared<OperatorNode>("BR(B>X_sll)_high_q2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return BR_B_Xsll(q2_high_bound.first, q2_high_bound.second, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16]); });
    br_high_tau->addChildren({m_b_1S, n_pref, m_tau_hat, alpha_s_mu_b, z, L_l_tau, L_b_5GeV, n_pref_dB0_ds, n_pref_delta_mb2, n_pref_delta_mb3, n_pref_delta_mc2, n_pref_delta_brems, n_pref_delta_em, L_b, m_c_hat, inv_alpha_em, m_D_hat, C9_eff_lookup, wilson_cache, n_delta_brems_B_lookup, n_f_17_lookup, n_f_19_lookup, n_f_27_lookup, n_f_29_lookup});

    roots.emplace(ObservableMapper::to_id(Observables::BR_B__Xs_mu_mu__LOW_Q2), br_low_mu);
    roots.emplace(ObservableMapper::to_id(Observables::BR_B__Xs_mu_mu__HIGH_Q2), br_high_mu);
    roots.emplace(ObservableMapper::to_id(Observables::BR_B__Xs_tau_tau__HIGH_Q2), br_high_tau);
}
