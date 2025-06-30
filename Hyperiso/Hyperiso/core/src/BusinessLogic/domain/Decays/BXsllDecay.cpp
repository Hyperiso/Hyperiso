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

scalar_t BXsllDecay::R_cc_cont(scalar_t s) {
    double s_re = s.real();
    if (s_re < 0 || s_re > 1) {
        LOG_WARN("In BXsllDecay::R_cc_cont(scalar_t s): s souldn't be outside [0, 1].");
    }

    if (s_re < .6) {
        return 0;
    } else if (s_re < 0.69) {
        -6.8 + 11.33 * s;
    } else {
        return 1.02;
    }
}

scalar_t BXsllDecay::breit_wigner(scalar_t s,
                                  scalar_t m_V,
                                  scalar_t br,
                                  scalar_t gamma_tot,
                                  scalar_t gamma_had)
{
    return br * gamma_tot * gamma_had / (pow(s - pow(m_V, 2), 2) + pow(m_V * gamma_tot, 2));
}

scalar_t BXsllDecay::R_cc(scalar_t s, scalar_t inv_alpha_em) {
    scalar_t R_cc_res;

    for (size_t k = 0; k < cc_res_mass.size(); k++) {
        R_cc_res += breit_wigner(s, cc_res_mass[k], cc_res_br[k], cc_res_width_tot[k], cc_res_width_had[k]);
    }

    return 9 * s * pow(inv_alpha_em, 2) * R_cc_res + R_cc_cont(s);
}

scalar_t BXsllDecay::g_ld(scalar_t z, scalar_t s, scalar_t inv_alpha_em, scalar_t m_D_hat) {
    auto f = [this, s, inv_alpha_em] (double sp) {
        return R_cc(sp, inv_alpha_em) / (sp * (sp - s));
    };

    double epsilon = 1e-4;
    scalar_t I_1 = integrate(f, 4 * pow(m_D_hat, 2), s - epsilon, 1e-3);
    scalar_t I_2 = integrate(f, s + epsilon, 1, 1e-3);

    return g(z, 0) + (s * (I_1 + I_2) + I * PI * R_cc(s, inv_alpha_em)) / 3.;
}

scalar_t BXsllDecay::C9_eff(scalar_t s, scalar_t L_mu, scalar_t mc_hat, scalar_t inv_alpha_em, scalar_t m_D_hat, QCDOrder order) {
    auto C = w_proxy->getAFR(WGroup::B, order);

    scalar_t g_0 = g(0, s);
    scalar_t g_1 = g(1, s);
    scalar_t g_mc = g_ld(mc_hat, s, inv_alpha_em, m_D_hat);

    return C[WCoef::C9]
	        +(-32./27.*C[WCoef::C1]-8./9.*C[WCoef::C2]-16./9.*C[WCoef::C3]+32./27.*C[WCoef::C4]-112./9.*C[WCoef::C5]+512./27.*C[WCoef::C6]) * L_mu
	        +4./3.*C[WCoef::C3]+64./9.*C[WCoef::C5]+64./27.*C[WCoef::C6]
	        +g_mc*(4./3.*C[WCoef::C1]+C[WCoef::C2]+6.*C[WCoef::C3]+60.*C[WCoef::C5])
	        +g_1*(-7./2.*C[WCoef::C3]-2./3.*C[WCoef::C4]-38.*C[WCoef::C5]-32./3.*C[WCoef::C6])
	        +g_0*(-1./2.*C[WCoef::C3]-2./3.*C[WCoef::C4]-8.*C[WCoef::C5]-32./3.*C[WCoef::C6]);
}

scalar_t BXsllDecay::F_17(scalar_t L_b, scalar_t z, scalar_t s) {
    scalar_t Lc = log(z) / 2;
	scalar_t Ls = log(s);
	scalar_t Lsb = log(1.-s);
	scalar_t Li2s = CLi2(s);
	scalar_t Li3s = CLi3(s);
	scalar_t Li3sb = CLi3(1.-s);
	scalar_t Li4s = CLi4(s);
	scalar_t sqrts = sqrt(s);
	scalar_t sqrt4s = sqrt(4-s);
	scalar_t as = asin(sqrts/2);
	scalar_t cl2 = Cl2(2*as);
	scalar_t cl3 = Cl3(2*as);

    std::array<scalar_t, 8> s_inv_pow;
    for (size_t k = 0; k < 8; k++) s_inv_pow[k] = pow(s, -(k + 1));

    std::array<scalar_t, 21> s_pow;
    for (size_t k = 0; k < 21; k++) s_pow[k] = pow(s, k + 1);

	scalar_t lz_sq = pow(lz,2);
	scalar_t lz_cube = lz*lz_sq;
	scalar_t lz_four = lz*lz_cube;
	scalar_t z_sqrt = pow(z,0.5);   
	scalar_t z_m1 = pow(z,-1.);      
	scalar_t z_m2 = z_m1*z_m1;      
	scalar_t z_sq = pow(z,2.);
	scalar_t z_cube = z*z_sq;
	scalar_t ms_p1_pow_2 = pow(1. - s,2.);
	scalar_t ms_p1_pow_3 = ms_p1_pow_2*(1.-s);
	scalar_t ms_p1_pow_4 = ms_p1_pow_2*ms_p1_pow_2*(1.-s);
	scalar_t s_m1_pow_m1 = pow(-1. + s,-1.);        
	scalar_t s_m4_s_sqrt = pow(-((-4 + s)*s),0.5);
	scalar_t s_m1_pow_m2 = s_m1_pow_m1*s_m1_pow_m1;
	scalar_t s_m1_pow_m3 = s_m1_pow_m1*s_m1_pow_m2;
	scalar_t s_m1_pow_m4 = s_m1_pow_m1*s_m1_pow_m3;
	scalar_t s_m1_pow_m5 = s_m1_pow_m1*s_m1_pow_m4;
	scalar_t s_m1_pow_m6 = s_m1_pow_m1*s_m1_pow_m5;
	scalar_t s_m1_pow_m7 = s_m1_pow_m1*s_m1_pow_m6;

	scalar_t res=0;

	if(abs(s) < 0.4) {
        for (kappa_coef k : kappa_17_low) {
            res += k.value * pow(s, k.i) * pow(Ls, k.j) * pow(z, k.l / 2.) * pow(Lc, k.m);
        }
        res += - 0.8559670781893005 * L_b;
    } else if(abs(s)<.900001) {
        res += -0.8559670781893004*L_b - 0.2962962962962963*s*pow(as,2)*s_m1_pow_m4 - 0.00823045267489712*Ls*(29. - 18.*Lsb*(-1. + s) - 47.*s)*s*
                s_m1_pow_m2 + 0.14814814814814814*Li2s*s*s_m1_pow_m1 - 0.07407407407407407*pow(Ls,2)*s_m1_pow_m3*s_pow[2] - 
                0.01646090534979424*as*sqrt4s*sqrts*s_m1_pow_m3*s_inv_pow[0]*(-4. + 9.*s - 15.*s_pow[1] + 4.*s_pow[2]) - 0.0013717421124828531*
                s_m1_pow_m4*(785. - 2.*s*(1585. + 12.*PI2) + 6.*(803. + 9.*PI2)*s_pow[1] - 2.*(1633. + 27.*PI2)*s_pow[2] + (833. + 18.*PI2)*s_pow[3]);

        res += z * (0.04938271604938271*s*pow(Ls,4)*s_m1_pow_m3 - 0.2962962962962963*Li2s*(-1. + 4.*s)* s_m1_pow_m2 + 0.09876543209876543*(-3. + s)*
                pow(Ls,3)*s_m1_pow_m2 - 1.*lz*(-0.2962962962962963*Ls*(2. + 2.*Lsb*(-1. + s) + s)*s_m1_pow_m2 - 0.5925925925925926*Li2s*s_m1_pow_m1 + 
                0.09876543209876543*(9. + PI2)*s_m1_pow_m1 + 0.2962962962962963*pow(Ls,2)*s_m1_pow_m1) - 0.35616500834358344*s_m1_pow_m3*
                (-5. + 5.*s + 2.*s_pow[1]) - 0.2962962962962963*Li3s*s_m1_pow_m3*(-1. - 5.*s + 4.*s_pow[1]) - 1.*pow(Ls,2)*(-0.5925925925925926*Li2s*s*
                s_m1_pow_m3 - 0.04938271604938271*s_m1_pow_m3*(12. + 3.*Lsb*(-5. + 3.*s) + s*(-15. + 2.*PI2) + 6.*s_pow[1])) - 
                1.*Ls*(-2.8493200667486676*s*s_m1_pow_m3 + 2.3703703703703702*Li3s*s*s_m1_pow_m3 + 0.2962962962962963*Li2s*s_m1_pow_m3*
                (3. + s - 2.*s_pow[1]) + 0.04938271604938271*s_m1_pow_m3*(12. + 21.*s + (-5. + 3.*s)*PI2 - 33.*s_pow[1] + 6.*Lsb*(1. - 5.*s + 4.*s_pow[1]))) - 
                1.7777777777777777*pow(as,2)*s_m1_pow_m4*(-1. - 3.*s_pow[1] + s_pow[2]) + 0.009876543209876543*s_m1_pow_m4*(-4.*(-1. + s)*s*pow(PI,4.) - 
                15.*(-1. + s)*(10. - (17. + 24.*Li4s)*s + 7.*s_pow[1]) + 5.*PI2*(-6. + 16.*s - 20.*s_pow[1] + 7.*s_pow[2])) - 
                0.8888888888888888*as*(1. + s)*s_m1_pow_m3*s_m4_s_sqrt);

        res += 0.19753086419753085*(2. + s)*PI2*s_inv_pow[0]*pow(z,1.5);

        res += (-0.2962962962962963*Li3s*(5. + 3.*s)*s_m1_pow_m3 + 0.2962962962962963*Li3sb*(1. + 6.*s)*s_m1_pow_m3 + 0.14814814814814814*pow(Ls,3)*
                s_m1_pow_m3 - 0.8888888888888888*as*(-3. + s)*sqrt4s*sqrts*s_m1_pow_m3*s_pow[1] + 0.35616500834358344*s_m1_pow_m4*
                (-5. + 8.*s + 3.*s_pow[1]) + 0.2962962962962963*Li2s*s_m1_pow_m3*s_inv_pow[0]*(3. - 10.*s + 8.*s_pow[1]) - 
                0.8888888888888888*pow(as,2)*s_m1_pow_m4*s_pow[1]*(2. + 9.*s - 6.*s_pow[1] + s_pow[2]) - 0.07407407407407407*pow(Ls,2)*s_m1_pow_m4*
                s_inv_pow[1]*(-4. + 7.*s + (5. + 6.*Lsb)*s_pow[1] - 1.*(25. + 8.*Lsb)*s_pow[2] + (11. + 2.*Lsb)*s_pow[3]) - 
                1.*lz_sq*(0.14814814814814814*Ls*(1. + 6.*s)*s_m1_pow_m3 - 0.07407407407407407*s_m1_pow_m2*s_inv_pow[1]*
                (1. - 2.*s + 8.*s_pow[1] + 10.*s_pow[2] - 4.*s_pow[3] + s_pow[4])) - 1.*lz*(-10.666666666666666*s*pow(as,2)*s_m1_pow_m4 - 
                0.2962962962962963*Li2s*(1. + 6.*s)*s_m1_pow_m3 - 0.14814814814814814*pow(Ls,2)*s_m1_pow_m2 - 
                0.14814814814814814*Ls*s_m1_pow_m3*s_inv_pow[1]*(2. - 3.*s + 2.*(1. + Lsb)*s_pow[1] + 3.*(-1. + 4.*Lsb)*s_pow[2]) + 
                0.024691358024691357*s_m1_pow_m4*s_inv_pow[1]*(9. - 42.*s + (81. - 2.*PI2)*s_pow[1] + 2.*(-6. + PI2)*s_pow[2] + 6.*(-23. + 2.*PI2)*s_pow[3]
                + 153.*s_pow[4] - 60.*s_pow[5] + 9.*s_pow[6])- 0.5925925925925926*as*s_m1_pow_m3*(-4. - 3.*s_pow[1] + s_pow[2])* s_m4_s_sqrt) - 
                1.*Ls*(-5.333333333333333*s*pow(as,2)*s_m1_pow_m4 - 0.2962962962962963*Li2s*(4. + s)*s_m1_pow_m3 - 
                0.024691358024691357*s_m1_pow_m3*s_inv_pow[1]*(-18. + (57. + 36.*Lsb)*s - 6.*(5. + 20.*Lsb + PI2)*s_pow[1] + (-93. + 96.*Lsb + 
                50.*PI2)*s_pow[2]) - 0.2962962962962963*as*s_m1_pow_m3*(-4. - 3.*s_pow[1] + s_pow[2])*s_m4_s_sqrt) + 
                0.012345679012345678*s_m1_pow_m4*s_inv_pow[1]*(2.*PI2*(-6. + 45.*s - 133.*s_pow[1] + 157.*s_pow[2] - 
                59.*s_pow[3] + 12.*s_pow[4] - 12.*s_pow[5] + 2.*s_pow[6]) + 3.*(7. - 78.*s + 7.*s_pow[6] + s_pow[4]*(155. - 32.*cl2*s_m4_s_sqrt) + 
                8.*s_pow[5]*(-7. + cl2*s_m4_s_sqrt) - 8.*s_pow[2]*(27. + 18.*cl3 + 4.*cl2*s_m4_s_sqrt) + s_pow[3]*(-62. + 24.*cl2*s_m4_s_sqrt) + 
                s_pow[1]*(243. + 32.*cl2*s_m4_s_sqrt))))*z_sq;

        res += 0.02633744855967078*PI2*s_inv_pow[2] * (3. + 14.*s + 3.*s_pow[1])*pow(z,2.5);

        res += (1.1851851851851851*Li3sb*s*s_m1_pow_m4 - 0.5925925925925926*Li3s*s_m1_pow_m3*s_inv_pow[0] + 0.19753086419753085*pow(Ls,3)*s_m1_pow_m3*
                s_inv_pow[0] - 0.7123300166871669*s_m1_pow_m4*s_inv_pow[0]*(1. - s + 2.*s_pow[1]) - 0.2962962962962963*Li2s*s_m1_pow_m4*s_inv_pow[1]*
                (4. - 14.*s + 18.*s_pow[1] - 11.*s_pow[2] + s_pow[3]) - 0.19753086419753085*cl2*s*sqrt4s*sqrts*s_m1_pow_m3*
                (-9. - s + 9.*s_pow[1] - 6.*s_pow[2] + s_pow[3]) - 0.06584362139917696*as*s*sqrt4s*sqrts*s_m1_pow_m3*(27. - s + 9.*s_pow[1] - 6.*s_pow[2] + 
                s_pow[3]) + 0.5925925925925926*s*pow(as,2)*s_m1_pow_m4*(-6. + 27.*s - 30.*s_pow[2] + 27.*s_pow[3] - 9.*s_pow[4] + s_pow[5]) - 
                1.*Ls*(3.5555555555555554*s*pow(as,2)*s_m1_pow_m4 - 0.5925925925925926*Li2s*s_m1_pow_m3*s_inv_pow[0] + 0.19753086419753085*as*s*sqrt4s*
                sqrts*s_m1_pow_m3*(-9. - s + 9.*s_pow[1] - 6.*s_pow[2] + s_pow[3]) + 0.00823045267489712*s_m1_pow_m4*s_inv_pow[3]*(-10. + 10.*s + 
                6.*(23. + 24.*Lsb)*s_pow[1] - 2.*(67. + 252.*Lsb + 30.*PI2)*s_pow[2] + (-163. + 648.*Lsb + 60.*PI2)*s_pow[3] - 12.*(-5. + 33.*Lsb + 
                6.*PI2)*s_pow[4] + 9.*(-5. + 4.*Lsb)*s_pow[5])) - 0.012345679012345678*pow(Ls,2)*s_m1_pow_m5*s_inv_pow[3]*(8. + 8.*s - 94.*s_pow[1] + 
                (94. - 24.*Lsb)*s_pow[2] + (52. + 48.*Lsb)*s_pow[3] - 8.*(19. + 3.*Lsb)*s_pow[4] + 72.*s_pow[5] - 12.*s_pow[6]) - 
                1.*lz_sq*(0.5925925925925926*Ls*s*s_m1_pow_m4 + 0.024691358024691357*s_m1_pow_m3*s_inv_pow[3]*(1. + 3.*s - 18.*s_pow[1] + 25.*s_pow[2] - 
                36.*s_pow[3] + 26.*s_pow[4] - 27.*s_pow[5] - 22.*s_pow[6] + 38.*s_pow[7] - 16.*s_pow[8] + 2.*s_pow[9])) - 
                1.*lz*(-1.1851851851851851*Li2s*s*s_m1_pow_m4 + 7.111111111111111*s*pow(as,2)*s_m1_pow_m4 + 0.2962962962962963*pow(Ls,2)*s_m1_pow_m3*
                s_inv_pow[0] + 0.3950617283950617*as*s*sqrt4s*sqrts*s_m1_pow_m3*(-9. - s + 9.*s_pow[1] - 6.*s_pow[2] + s_pow[3]) + 0.04938271604938271*Ls*
                s_m1_pow_m4*s_inv_pow[3]*(2. + 4.*s - 27.*s_pow[1] + 22.*s_pow[2] + 14.*s_pow[3] - 12.*(1. + 2.*Lsb)*s_pow[4] + 9.*s_pow[5]) + 
                0.00823045267489712*s_m1_pow_m3*s_inv_pow[3]*(-5. + 57.*s_pow[1] - 46.*s_pow[2] - 85.*s_pow[3] + 164.*s_pow[4] - 73.*s_pow[5] + 
                102.*s_pow[6] - 40.*s_pow[7] - 4.*s_pow[8] + 2.*s_pow[9])) - 0.0013717421124828531*s_m1_pow_m4*s_inv_pow[3]*
                (-19. + 46.*s + 1365.*s_pow[1] - 4435.*s_pow[2] + 3757.*s_pow[3] + (610. - 2592.*cl3)*s_pow[4] - 3103.*s_pow[5] + 2113.*s_pow[6] + 
                594.*s_pow[7] - 1548.*s_pow[8] + 714.*s_pow[9] - 94.*s_pow[10] + 6.*PI2*(6. + 12.*s - 195.*s_pow[1] + 442.*s_pow[2] - 
                434.*s_pow[3] + 124.*s_pow[4] + 137.*s_pow[5] - 16.*s_pow[6] - 120.*s_pow[7] + 108.*s_pow[8] - 36.*s_pow[9] + 4.*s_pow[10])))*z_cube;

        res += 0.002257495590828924*PI2*s_inv_pow[4]* (15. + 108.*s + 314.*s_pow[1] + 108.*s_pow[2] + 15.*s_pow[3])*pow(z,3.5);

        res += (1.1872166944786116*s_m1_pow_m3*s_inv_pow[1] - 0.9876543209876543*Li3s*s_m1_pow_m3*s_inv_pow[1] + 0.3292181069958848*pow(Ls,3)*s_m1_pow_m3*
                s_inv_pow[1] + 0.2962962962962963*cl2*s*sqrt4s*sqrts*(8. + 6.*s - 6.*s_pow[1] + s_pow[2]) - 
                0.04938271604938271*Li2s*s_m1_pow_m5*s_inv_pow[2] *(-50. + 230.*s - 420.*s_pow[1] + 380.*s_pow[2] - 146.*s_pow[3] + 15.*s_pow[4]) + 
                0.04938271604938271*as*s*sqrt4s*sqrts*s_m1_pow_m3*(16. + 72.*s + 18.*s_pow[1] - 203.*s_pow[2] + 189.*s_pow[3] - 63.*s_pow[4] + 
                7.*s_pow[5]) - 0.8888888888888888*pow(as,2)*s_m1_pow_m4*s_pow[1]*(36. - 60.*s + 99.*s_pow[2] - 112.*s_pow[3] + 
                54.*s_pow[4] - 12.*s_pow[5] + s_pow[6]) - 1.*Ls*(-0.9876543209876543*Li2s*s_m1_pow_m3*s_inv_pow[1] - 0.2962962962962963*as*s*sqrt4s*
                sqrts*(8. + 6.*s - 6.*s_pow[1] + s_pow[2]) + 0.0013717421124828531*s_m1_pow_m5*s_inv_pow[5]*(21. - 10.*s - 307.*s_pow[1] - 
                1.*(283. + 1800.*Lsb)*s_pow[2] + 5.*(437. + 1656.*Lsb + 120.*PI2)*s_pow[3] - 2.*(1081. + 7560.*Lsb + 600.*PI2)*s_pow[4] + 
                (253. + 13680.*Lsb + 600.*PI2)*s_pow[5] + (861. - 5256.*Lsb)*s_pow[6] + 6.*(49. + 90.*Lsb)*s_pow[7] - 78.*s_pow[8]))
                + 0.00823045267489712*pow(Ls,2)*s_m1_pow_m6*s_inv_pow[5]*(6. + 14.*s - 30.*s_pow[1] - 255.*s_pow[2] + (729. - 60.*Lsb)*s_pow[3] + 
                45.*(-13. + 4.*Lsb)*s_pow[4] - 15.*(11. + 12.*Lsb)*s_pow[5] + (541. + 60.*Lsb)*s_pow[6] - 270.*s_pow[7] + 45.*s_pow[8]) + 
                0.00411522633744856*lz_sq*s_inv_pow[5]*(3. + 25.*s + 90.*s_pow[1] - 60.*s_pow[2] + 5.*s_pow[3] - 3.*s_pow[4] - 144.*s_pow[5] - 
                252.*s_pow[6] + 72.*s_pow[7] + 288.*s_pow[8] - 144.*s_pow[9] + 18.*s_pow[10]) - 1.*lz*(0.49382716049382713*pow(Ls,2)*s_m1_pow_m3*
                s_inv_pow[1] - 0.5925925925925926*as*s*sqrt4s*sqrts*(8. + 6.*s - 6.*s_pow[1] + s_pow[2]) + 0.01646090534979424*Ls*s_m1_pow_m5*s_inv_pow[5]*
                (-3. - 10.*s + 5.*s_pow[1] + 185.*s_pow[2] - 425.*s_pow[3] + 343.*s_pow[4] - 125.*s_pow[5] + 3.*s_pow[6] - 15.*s_pow[7] + 15.*s_pow[8]) + 
                0.00013717421124828533*s_m1_pow_m4*s_inv_pow[5]*(105. + 55.*s - 1480.*s_pow[1] + 90.*s_pow[2] + 10470.*s_pow[3] - 25782.*s_pow[4] + 
                43242.*s_pow[5] - 52678.*s_pow[6] + 52987.*s_pow[7] - 14283.*s_pow[8] - 23976.*s_pow[9] + 3240.*s_pow[10] + 
                28080.*s_pow[11] - 22680.*s_pow[12] + 6480.*s_pow[13] - 630.*s_pow[14])) + 2.2862368541380887e-6*s_m1_pow_m5*s_inv_pow[5]*
                (-2775. + 5550.*s + 34425.*s_pow[1] + 1.48905e6*s_pow[2] - 7.4002e6*s_pow[3] + 1.4083748e7*s_pow[4] - 1.4200216e7*s_pow[5] + 
                9.34426e6*s_pow[6] - 5.568515e6*s_pow[7] + 141850.*s_pow[8] + 4.948647e6*s_pow[9] - 27324.*s_pow[10] - 8.9613e6*s_pow[11] + 
                9.8388e6*s_pow[12] - 4.6764e6*s_pow[13] + 1.03725e6*s_pow[14] - 86850.*s_pow[15] + 1200.*PI2*(9. + 30.*s - 15.*s_pow[1] - 1245.*s_pow[2] + 
                4325.*s_pow[3] - 6463.*s_pow[4] + 5357.*s_pow[5] - 2657.*s_pow[6] + 295.*s_pow[7] + 1405.*s_pow[8] - 1014.*s_pow[9] - 1782.*s_pow[10] + 
                3798.*s_pow[11] - 2988.*s_pow[12] + 1188.*s_pow[13] - 234.*s_pow[14] + 18.*s_pow[15])))*pow(z,4.);

        res += 0.0005374989501973629*PI2*s_inv_pow[6]*(35. + 330.*s + 1389.*s_pow[1] + 3212.*s_pow[2] + 1389.*s_pow[3] + 330.*s_pow[4] + 35.*s_pow[5])*pow(z,4.5);

        res += (2.493155058405084*s_m1_pow_m3*s_inv_pow[2] - 2.074074074074074*Li3s*s_m1_pow_m3*s_inv_pow[2] + 0.691358024691358*pow(Ls,3)*s_m1_pow_m3*
                s_inv_pow[2]  - 0.5925925925925926*cl2*sqrt4s*sqrts*s_pow[1]*(-20. - 11.*s + 24.*s_pow[1] - 9.*s_pow[2] + s_pow[3]) - 
                0.019753086419753086*as*sqrt4s*sqrts*s_m1_pow_m3*s_pow[1]*(630. - 1853.*s - 99.*s_pow[1] + 4458.*s_pow[2] - 5217.*s_pow[3] + 2538.*s_pow[4] - 
                564.*s_pow[5] + 47.*s_pow[6]) + 0.09876543209876543*Li2s*s_m1_pow_m6*s_inv_pow[3]*(-63. + 355.*s - 826.*s_pow[1] + 1008.*s_pow[2] - 
                665.*s_pow[3] + 213.*s_pow[4] - 7.*s_pow[6] + s_pow[7]) + 1.7777777777777777*pow(as,2)*s_m1_pow_m4*s_pow[1]*
                (10. - 90.*s + 165.*s_pow[1] - 333.*s_pow[3] + 445.*s_pow[4] - 275.*s_pow[5] + 90.*s_pow[6] - 15.*s_pow[7] + s_pow[8]) - 
                1.*Ls*(-2.074074074074074*Li2s*s_m1_pow_m3*s_inv_pow[2]  + 0.5925925925925926*as*sqrt4s*sqrts*s_pow[1]*(-20. - 11.*s + 24.*s_pow[1] - 
                9.*s_pow[2] + s_pow[3]) + 0.0000823045267489712*s_m1_pow_m6*s_inv_pow[7]*(-162. + 36.*s + 2914.*s_pow[1] + 5862.*s_pow[2] + 
                30.*(-767. + 2520.*Lsb)*s_pow[3] - 8.*(2131. + 53250.*Lsb + 2625.*PI2)*s_pow[4] + 8.*(12113. + 123900.*Lsb + 7875.*PI2)*s_pow[5] - 
                18.*(5693. + 67200.*Lsb + 3500.*PI2)*s_pow[6] + 3.*(15101. + 266000.*Lsb + 7000.*PI2)*s_pow[7] - 50.*(259. + 5112.*Lsb)*s_pow[8]
                + 9600.*s_pow[9] + 150.*(11. + 56.*Lsb)*s_pow[10] - 25.*(25. + 48.*Lsb)*s_pow[11])) - 0.0024691358024691358*pow(Ls,2)*s_m1_pow_m7*
                s_inv_pow[7]*(12. + 42.*s - 28.*s_pow[1] - 308.*s_pow[2] - 1533.*s_pow[3] + (8318. - 420.*Lsb)*s_pow[4] + 6.*(-2237. + 280.*Lsb)*s_pow[5] - 
                56.*(-143. + 45.*Lsb)*s_pow[6] + 12.*(149. + 140.*Lsb)*s_pow[7] - 1.*(4507. + 420.*Lsb)*s_pow[8] + 1540.*s_pow[9] + 140.*s_pow[10] - 
                160.*s_pow[11] + 20.*s_pow[12]) - 0.0012345679012345679*lz_sq*s_inv_pow[7]*(-6. - 63.*s - 301.*s_pow[1] - 840.*s_pow[2] + 
                630.*s_pow[3] - 99.*s_pow[4] + 153.*s_pow[5] + 6.*s_pow[6] + 1200.*s_pow[7] + 1440.*s_pow[8] + 3360.*s_pow[9] - 2400.*s_pow[10] - 
                5400.*s_pow[11] + 4800.*s_pow[12] - 1320.*s_pow[13] + 120.*s_pow[14]) - 1.*lz*(1.037037037037037*pow(Ls,2)*s_m1_pow_m3*
                s_inv_pow[2]  + 1.1851851851851851*as*sqrt4s*sqrts*s_pow[1]*(-20. - 11.*s + 24.*s_pow[1] - 9.*s_pow[2] + s_pow[3]) + 
                0.0049382716049382715*Ls*s_m1_pow_m6*s_inv_pow[7]*(6. + 27.*s + 13.*s_pow[1] - 141.*s_pow[2] - 1380.*s_pow[3] + 
                5549.*s_pow[4] - 8207.*s_pow[5] + 6207.*s_pow[6] - 2809.*s_pow[7] + 1080.*s_pow[8] - 30.*s_pow[9] - 30.*s_pow[10] + 35.*s_pow[11]) + 
                0.000011757789535567314*s_m1_pow_m5*s_inv_pow[7]*(-567. - 441.*s + 9758.*s_pow[1] + 30275.*s_pow[2] - 166495.*s_pow[3] + 
                20923.*s_pow[4] + 866054.*s_pow[5] - 1.788427e6*s_pow[6] + 1.255755e6*s_pow[7] + 458890.*s_pow[8] - 1.258076e6*s_pow[9] + 
                55832.*s_pow[10] - 524481.*s_pow[11] + 2.826e6*s_pow[12] - 901320.*s_pow[13] - 5.25168e6*s_pow[14] + 7.8771e6*s_pow[15] - 
                5.0799e6*s_pow[16] + 1.7073e6*s_pow[17] - 290640.*s_pow[18] + 19740.*s_pow[19])) - 1.3997368494722993e-8*s_m1_pow_m6*s_inv_pow[7]*
                (-161406. + 318843.*s + 1.910657e6*s_pow[1] - 738969.*s_pow[2] + 6.33135615e8*s_pow[3] - 3.654733992e9*s_pow[4] + 8.580730746e9*s_pow[5] - 
                1.0387902666e10*s_pow[6] + 6.272402572e9*s_pow[7] - 4.35513965e8*s_pow[8] - 2.310469467e9*s_pow[9] + 2.213537811e9*s_pow[10] - 
                4.191847021e9*s_pow[11] + 8.659496442e9*s_pow[12] - 5.00515056e9*s_pow[13] - 1.057709688e10*s_pow[14] + 2.294859924e10*s_pow[15] - 
                2.03364504e10*s_pow[16] + 1.00250472e10*s_pow[17] - 2.83732932e9*s_pow[18] + 4.2896364e8*s_pow[19] - 2.674812e7*s_pow[20] + 58800.*PI2*
                (18. + 81.*s + 39.*s_pow[1] - 423.*s_pow[2] - 9810.*s_pow[3] + 47729.*s_pow[4] - 95467.*s_pow[5] + 103887.*s_pow[6] - 
                63289.*s_pow[7] + 15730.*s_pow[8] + 8388.*s_pow[9] - 21404.*s_pow[10] + 46029.*s_pow[11] - 50988.*s_pow[12] - 
                20160.*s_pow[13] + 133320.*s_pow[14] - 179760.*s_pow[15] + 130200.*s_pow[16] -  56400.*s_pow[17] + 14520.*s_pow[18] - 
                2040.*s_pow[19] + 120.*s_pow[20])))*pow(z,5.);
    } else {
        for (kappa_coef k : kappa_17_high) {
            res += k.value * pow(s, k.i) * pow(Ls, k.j) * pow(z, k.l / 2.) * pow(Lc, k.m);
        }
        res += -0.8559670781893005 * L_b;
    }

	return res;
}
