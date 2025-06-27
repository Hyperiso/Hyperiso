#include "BXsllDecay.h"

scalar_t BXsllDecay::alpha_s(scalar_t mu) {
    return ObsQCDProxy()(AlphasConfig{mu, MassType::POLE, MassType::POLE});
}

scalar_t BXsllDecay::f(scalar_t z) {
    return 1 - 8 * z + 8 * pow(z, 3) - pow(z, 4) - 12 * pow(z, 2) * log(z);
}

scalar_t BXsllDecay::h(scalar_t z) {
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

scalar_t BXsllDecay::kappa(scalar_t f, scalar_t h, scalar_t alpha_s_mu_b) {
    return 1 - 2 * alpha_s_mu_b * h / (3 * PI * f);
}

scalar_t BXsllDecay::m_hat(scalar_t m, scalar_t mb) {
    return m / mb;
}

scalar_t BXsllDecay::z(scalar_t mc_hat) {
    return pow(mc_hat, 2);
}

scalar_t BXsllDecay::f_7(scalar_t s) {
    return 1./6./(s-1.)/(s-1.)*(24.*(1.+13.*s-4.*s*s)*Li2(sqrt(s))+12.*(1.-17.*s+6.*s*s)*Li2(s)+6.*s*(6.-7.*s)*log(s)
	+24.*(1.-s)*(1.-s)*log(s)*log(1.-s)+12.*(-13.+16.*s-3.*s*s)*(log(1.-sqrt(s))-log(1.-s))
	+39.-2.*PI2+252.*s-26.*PI2*s+21.*s*s+8.*PI2*s*s-180.*sqrt(s)-132.*s*sqrt(s));
}

scalar_t BXsllDecay::f_9(scalar_t s) {
    return -1./6./(s-1.)/(s-1.)*(48.*s*(-5.+2.*s)*Li2(sqrt(s))+24.*(-1.+7.*s-3.*s*s)*Li2(s)+6.*s*(-6.+7.*s)*log(s)
	-24.*(1.-s)*(1.-s)*log(s)*log(1.-s)+24.*(5.-7.*s+2.*s*s)*(log(1.-sqrt(s))-log(1.-s))
	-21.-156.*s+20.*PI2*s+9.*s*s-8.*PI2*s*s+120.*sqrt(s)+48.*s*sqrt(s));
}

scalar_t BXsllDecay::tau_77(scalar_t s) {
    return -2./9./(2.+s)*(2.*(1.-s)*(1.-s)*log(1.-s)+6.*s*(2.-2.*s-s*s)/(1.-s)/(1.-s)*log(s)+(11.-7.*s-10.*s*s)/(1.-s));
}

scalar_t BXsllDecay::tau_99(scalar_t s) {
    return -4./9./(1.+2.*s)*(2.*(1.-s)*(1.-s)*log(1.-s)+3.*s*(1.+s)*(1.-2.*s)/(1.-s)/(1.-s)*log(s)+3.*(1.-3.*s*s)/(1.-s));
}

scalar_t BXsllDecay::tau_79(scalar_t s) {
    return -4.*(1.-s)*(1.-s)/9./s*log(1.-s)-4.*s*(3.-2.*s)*log(s)/9./(1.-s)/(1.-s)-2./9.*(5.-3.*s)/(1.-s);
}

scalar_t BXsllDecay::tau_710(scalar_t s) {
    return -5./2.+1./3./(1.-3.*s)-1./3.*s*(6.-7.*s)*log(s)/(1.-s)/(1.-s)-1./9.*(3.-7.*s+4.*s*s)*log(1.-s)/s+f_7(s)/3.;
}

scalar_t BXsllDecay::tau_910(scalar_t s) {
    return -5./2.+1./3./(1.-s)-1./3.*s*(6.-7.*s)*log(s)/(1.-s)/(1.-s)-2./9.*(3.-5.*s+2.*s*s)*log(1.-s)/s+f_9(s)/3.;
}

scalar_t BXsllDecay::sigma(scalar_t s) {
    return -4./3.*Li2(s)-2./3.*log(s)*log(1.-s)-2./9.*PI2-log(1.-s)-2./9.*(1.-s)*log(1.-s);
}

scalar_t BXsllDecay::sigma_9(scalar_t s) {
    return sigma(s)+1.5;
}

scalar_t BXsllDecay::sigma_7(scalar_t s, scalar_t L_mu) {
    return sigma(s)+1./6.-8./3.*L_mu;
}

scalar_t BXsllDecay::g(scalar_t z, scalar_t s) {
    scalar_t z2=z*z;

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

scalar_t BXsllDecay::C7_eff() {
    return this->w_proxy->getFR(WGroup::B, WCoef::C7, this->w_config.order);
}
