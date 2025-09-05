#include "BKsllDecay.h"

std::shared_ptr<OperatorNode> BKstarllDecay::load_FF_params() {
    auto ff_params = std::make_shared<OperatorNode>("FF_params", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 0; });
    int ff_id = (int)(this->cfg.ff_src) + 1;
    int sse_order = this->cfg.ff_src == BKstarllConfig::FF_Src::HLMW ? 1 : 2;
    
    auto get_m = [ff_id] (int i) { return ObsParameterProxy()(ParamId{ParameterType::DECAY, "B_Ksll", {ff_id, 0, i}}); };
    cache.m_R[FF::A0] = get_m(1);
    cache.m_R[FF::V] = cache.m_R[FF::T1] = get_m(2);
    cache.m_R[FF::A1] = cache.m_R[FF::A12] = cache.m_R[FF::T2] = cache.m_R[FF::T23] = get_m(3);

    for (int j = 1; j <= 3; j++) {
        ff_params->addChild(std::make_shared<ParameterNode>(ParamId{ParameterType::DECAY, "B_Ksll", {ff_id, 0, j}}));
    }

    for (int i = 1; i <= 7; i++) {
        for (int j = 0; j <= sse_order; j++) {
            ParamId PId {ParameterType::DECAY, "B_Ksll", {ff_id, i, j}};
            cache.alpha_ai[(FF)(i - 1)][j] = ObsParameterProxy()(PId);
            ff_params->addChild(std::make_shared<ParameterNode>(PId));
        }
    }

    FF_loaded = true;
    return ff_params;
}

double BKstarllDecay::z(double t) {
    double a = std::sqrt(cache.tp - t);
    double b = std::sqrt(cache.tp - cache.t0);
    return (a - b) / (a + b);
}

double BKstarllDecay::pole(double q2, double m_R) {
    return 1 / (1 - q2 / std::pow(m_R, 2));
}

complex_t BKstarllDecay::h(double q2, double m_q, double mu_b) {
    if(fpeq(m_q, 0.)) return 4./9.*(2./3.+I*PI+std::log(mu_b * mu_b / q2));
	
	double z=4.*m_q*m_q/q2;
    double L = 2 * std::log(m_q / mu_b);
	
	if(z>1.) return -4./9.*(L-2./3.-z)
	-4./9.*(2.+z)*sqrt(z-1.)*atan(1./sqrt(z-1.));
	
	else return -4./9.*(L-2./3.-z)
	-4./9.*(2.+z)*sqrt(1.-z)*(log((1.+sqrt(1.-z))/sqrt(z))-I*PI/2.);
}

double BKstarllDecay::phi_Kstar(double u, double a1, double a2) {
    double x=2.*u-1.;
	double C1=3.*x;
	double C2=-1.5+15./2.*x*x;

	return 6.*u*(1.-u)*(1.+a1*C1+a2*C2);
}

complex_t BKstarllDecay::B_0(double s, double m_q) {
    double epsilon=1.e-10;
	return -2.*std::sqrt(4.*(m_q*m_q-I*epsilon)/s-1.)*std::atan(1./std::sqrt(4.*(m_q*m_q-I*epsilon)/s-1.));
}

complex_t BKstarllDecay::L_1(complex_t x) {
    return std::log((x-1.)/x)*std::log(1.-x)-PI2/6.+CLi2(x/(x-1.));
}

complex_t BKstarllDecay::I_1(double u, double m_q, double q2, double m_Bd) {
    if(m_q==0.) return 1.;
	
	double epsilon=1.e-10;
    double mq2 = m_q * m_q;

	complex_t xp=0.5+std::sqrt(0.25-(mq2-I*epsilon)/((1.-u)*m_Bd*m_Bd+u*q2));
	complex_t xm=0.5-std::sqrt(0.25-(mq2-I*epsilon)/((1.-u)*m_Bd*m_Bd+u*q2));
	complex_t yp=0.5+std::sqrt(0.25-(mq2-I*epsilon)/q2);
	complex_t ym=0.5-std::sqrt(0.25-(mq2-I*epsilon)/q2);

	return 1.+2.*mq2/(1.-u)/(m_Bd*m_Bd-q2)*(L_1(xp)+L_1(xm)-L_1(yp)-L_1(ym));
}

complex_t BKstarllDecay::Y(double q2) {
    return h(q2, cache.m_c_pole, w_config.hadronic_scale) * (4./3. * cache.C[WCoef::C1] + cache.C[WCoef::C2] + 6. * cache.C[WCoef::C3] + 60. * cache.C[WCoef::C5])
            - 0.5 * h(q2, cache.m_b_pole, w_config.hadronic_scale) * (7. * cache.C[WCoef::C3] + 4./3. * cache.C[WCoef::C4] + 76. * cache.C[WCoef::C5] + 64./3. * cache.C[WCoef::C6])
            - 0.5 * h(q2, 0, w_config.hadronic_scale) * (cache.C[WCoef::C3] + 4./3. * cache.C[WCoef::C4] + 16. * cache.C[WCoef::C5] + 64./3. * cache.C[WCoef::C6])
            + 4./3. * cache.C[WCoef::C3] + 64./9. * cache.C[WCoef::C5] + 64./27. * cache.C[WCoef::C6];
}

complex_t BKstarllDecay::t_perp(double u, double m_q, double q2, double E_Kstar, double m_Bd) {
    if(fpeq(q2, 0.)) {
		if (fpeq(m_q, 0.)) return 4./(1.-u);
		double epsilon=1.e-10;
		complex_t xp=0.5+std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*m_Bd*m_Bd));
		complex_t xm=0.5-std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*m_Bd*m_Bd));
		return 4./(1.-u)*(1.+2.*m_q*m_q/(1.-u)/m_Bd/m_Bd*(L_1(xp)+L_1(xm)));
	}
	else return 2.*m_Bd/(1.-u)/E_Kstar*I_1(u,m_q,q2,m_Bd)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B_0((1.-u)*m_Bd*m_Bd+u*q2,m_q)-B_0(q2,m_q));
}

complex_t BKstarllDecay::t_par(double u, double m_q, double q2, double E_Kstar, double m_Bd) {
    return 2.*m_Bd/(1.-u)/E_Kstar*I_1(u,m_q,q2,m_Bd)+((1.-u)*m_Bd*m_Bd+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B_0((1.-u)*m_Bd*m_Bd+u*q2,m_q)-B_0(q2,m_q));
}

complex_t BKstarllDecay::F_27_u(double s_hat, double l) {
    double z=4./s_hat;

	complex_t A=
	208./243.*l+4.*s_hat/27./(1.-s_hat)*(Li2(s_hat)+log(s_hat)*log(1.-s_hat))
	+1./729./(1.-s_hat)/(1.-s_hat)*(6.*s_hat*(29.-47.*s_hat)*log(s_hat)+785.-1600.*s_hat+833.*s_hat*s_hat+6.*PI*I*(20.-49.*s_hat+47.*s_hat*s_hat))
	-2./243./pow(1.-s_hat,3.)*(2.*std::sqrt(z-1.)*(-4.+9.*s_hat-15.*s_hat*s_hat+4.*s_hat*s_hat*s_hat)*(PI/2.-std::atan(std::sqrt(z-1.)))+9.*s_hat*s_hat*s_hat*log(s_hat)*log(s_hat)+18.*PI*I*s_hat*(1.-2.*s_hat)*log(s_hat))
	+2.*s_hat/243./pow(1.-s_hat,4.)*(36.*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)+PI2*(-4.+9.*s_hat-9.*s_hat*s_hat+3.*s_hat*s_hat*s_hat));
	
	return -6.*A;
}

complex_t BKstarllDecay::F_19_u(double s_hat, double l) {
    double z=4./s_hat;
	
	complex_t x1=0.5+0.5*I*std::sqrt(z-1.);
	complex_t x2=0.5-0.5*I*std::sqrt(z-1.);
	complex_t x3=0.5+0.5*I/std::sqrt(z-1.);
	complex_t x4=0.5-0.5*I/std::sqrt(z-1.);

	complex_t B=
	8./243./s_hat*(-2.*(4.-34.*s_hat-17.*PI*I*s_hat)*l+8.*s_hat*pow(2.*l,2.)-17.*s_hat*log(s_hat)*2.*l)
	+(2.+s_hat)*std::sqrt(z-1.)/729./s_hat*(48.*2.*l*(PI/2.-std::atan(std::sqrt(z-1.)))-18.*PI*std::log(z-1.)+3.*I*std::log(z-1.)*std::log(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*PI2*I+6.*I*(-9.*std::log(x1)*std::log(x1)+std::log(x2)*std::log(x2)-2.*std::log(x4)*std::log(x4)+6.*std::log(x1)*std::log(x2)-4.*std::log(x1)*std::log(x3)+8.*std::log(x1)*std::log(x4))
	-12.*PI*(2.*std::log(x1)+std::log(x3)+std::log(x4)))
	-2./243./s_hat/(1.-s_hat)*(4.*s_hat*(-8.+17.*s_hat)*(Li2(s_hat)+log(s_hat)*log(1.-s_hat))
	+3.*(2.+s_hat)*(3.-s_hat)*std::log(x2/x1)*std::log(x2/x1)+12.*PI*(-6.-s_hat+s_hat*s_hat)*(PI/2.-std::atan(std::sqrt(z-1.))))
	+2./2187./s_hat/(1.-s_hat)/(1.-s_hat)*(-18.*s_hat*(120.-211.*s_hat+73.*s_hat*s_hat)*log(s_hat)-288.-8.*s_hat+934.*s_hat*s_hat-692.*s_hat*s_hat*s_hat+18.*PI*I*s_hat*(82.-173.*s_hat+73.*s_hat*s_hat))
	-4./243./s_hat/pow(1.-s_hat,3.)*(-2.*std::sqrt(z-1.)*(4.-3.*s_hat-18.*s_hat*s_hat+16.*s_hat*s_hat*s_hat-5.*pow(s_hat,4.))*(PI/2.-std::atan(std::sqrt(z-1.)))-9.*s_hat*s_hat*s_hat*log(s_hat)*log(s_hat)+2.*PI*I*s_hat*(8.-33.*s_hat+51.*s_hat*s_hat-17.*s_hat*s_hat*s_hat)*log(s_hat))
	+2./729./s_hat/pow(1.-s_hat,4.)*(72.*(3.-8.*s_hat+2.*s_hat*s_hat)*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)-PI2*(54.-53.*s_hat-286.*s_hat*s_hat+612.*pow(s_hat,3.)-446.*pow(s_hat,4.)+113.*pow(s_hat,5.)));
		
	complex_t C=-16./81.*(log(s_hat)-2.*l)+428./243.-64./27.*ZETA3+16./81.*PI*I;
		
	return B+4.*C;
}

complex_t BKstarllDecay::F_29_u(double s_hat, double l) {
    double z=4./s_hat;

	complex_t x1=0.5+0.5*I*std::sqrt(z-1.);
	complex_t x2=0.5-0.5*I*std::sqrt(z-1.);
	complex_t x3=0.5+0.5*I/std::sqrt(z-1.);
	complex_t x4=0.5-0.5*I/std::sqrt(z-1.);

	complex_t B=
	8./243./s_hat*(-2.*(4.-34.*s_hat-17.*PI*I*s_hat)*l+8.*s_hat*pow(2.*l,2.)-17.*s_hat*log(s_hat)*2.*l)
	+(2.+s_hat)*std::sqrt(z-1.)/729./s_hat*(48.*2.*l*(PI/2.-std::atan(std::sqrt(z-1.)))-18.*PI*std::log(z-1.)+3.*I*std::log(z-1.)*std::log(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*PI2*I+6.*I*(-9.*std::log(x1)*std::log(x1)+std::log(x2)*std::log(x2)-2.*std::log(x4)*std::log(x4)+6.*std::log(x1)*std::log(x2)-4.*std::log(x1)*std::log(x3)+8.*std::log(x1)*std::log(x4))
	-12.*PI*(2.*std::log(x1)+std::log(x3)+std::log(x4)))
	-2./243./s_hat/(1.-s_hat)*(4.*s_hat*(-8.+17.*s_hat)*(Li2(s_hat)+log(s_hat)*log(1.-s_hat))
	+3.*(2.+s_hat)*(3.-s_hat)*std::log(x2/x1)*std::log(x2/x1)+12.*PI*(-6.-s_hat+s_hat*s_hat)*(PI/2.-std::atan(std::sqrt(z-1.))))
	+2./2187./s_hat/(1.-s_hat)/(1.-s_hat)*(-18.*s_hat*(120.-211.*s_hat+73.*s_hat*s_hat)*log(s_hat)-288.-8.*s_hat+934.*s_hat*s_hat-692.*s_hat*s_hat*s_hat+18.*PI*I*s_hat*(82.-173.*s_hat+73.*s_hat*s_hat))
	-4./243./s_hat/pow(1.-s_hat,3.)*(-2.*std::sqrt(z-1.)*(4.-3.*s_hat-18.*s_hat*s_hat+16.*s_hat*s_hat*s_hat-5.*pow(s_hat,4.))*(PI/2.-std::atan(std::sqrt(z-1.)))-9.*s_hat*s_hat*s_hat*log(s_hat)*log(s_hat)+2.*PI*I*s_hat*(8.-33.*s_hat+51.*s_hat*s_hat-17.*s_hat*s_hat*s_hat)*log(s_hat))
	+2./729./s_hat/pow(1.-s_hat,4.)*(72.*(3.-8.*s_hat+2.*s_hat*s_hat)*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)-PI2*(54.-53.*s_hat-286.*s_hat*s_hat+612.*pow(s_hat,3.)-446.*pow(s_hat,4.)+113.*pow(s_hat,5.)));
		
	complex_t C=-16./81.*(log(s_hat)-2.*l)+428./243.-64./27.*ZETA3+16./81.*PI*I;	
	return -6.*B+3.*C;
}

complex_t BKstarllDecay::A_Seidel(double s, double m_b, double mu_b) {
    double shat = s/m_b/m_b;
    double z = 4./shat;

    // In the limit s -> 0 limit
    if (fabs(s) <= 1e-6) {
    return (1.1426611796982167 + 0.517134593183505 * I) +
    shat * ((-0.3221817635475234 - 0.23271056693257725 * I) +
    shat * ((-0.16999092571433885 + 0.23271056693257727 * I) +
    ((-0.13262865040070346 + 0.6981317007977318 * I) -
    (0.11612865336448869 - 1.1635528346628865 * I) * shat) * shat)) -
    0.8559670781893004 * log(m_b/mu_b) -
    shat * ((0.23868312757201646 - 0.4654211338651545 * I) +
    shat * ((-0.05761316872427984 - 0.46542113386515455 * I) +
    (-0.27983539094650206 - (0.4773662551440329 -
    0.9308422677303091 * I) * shat) * shat) +
    (-0.07407407407407407 - 0.22222222222222224 * shat) * (-shat * shat / 2));
    }

    complex_t mu_b_term = -104. / 243. * 2. * log(m_b/mu_b);
    if (fabs(shat - 1.0) < 1e-2) {
        const complex_t c0 = (997. + 18. * sqrt(3.) * PI) / 1458. + I * (64./243. * PI);
        const complex_t c1 = (215. + 9. * sqrt(3.) * PI) / 1215. + I * (-1./27. * PI);
        const complex_t c2 = (95. + 12. * sqrt(3.) * PI) / 2430. + I * (-7./405. * PI);

        return mu_b_term + c0 + c1 * (1.-shat) + c2 * pow(1.-shat, 2);
    }

    // TODO : Check consistency of std::log with gsl::Li2 (see old comment below)
    /* In the A expression, for (CLi2(shat)+log(shat)*std::log(1.-shat)), the real part is calculated correctly
     * but the imaginary part above the branch cut is not (in principle it should get cancelled).
     * I fixed the CLi2 to use the same branch cut as std::log, so it can be used directly */
    complex_t Li2log_term = CLi2(shat)+std::log(shat)*std::log(1.-shat);

    return -104./243.*log(m_b*m_b/mu_b/mu_b)+4.*shat/27./(1.-shat)*Li2log_term
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*std::log(shat)+785.-1600.*shat+833.*shat*shat+6.*PI*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*std::sqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(PI/2.-std::atan(std::sqrt(z-1.)))+9.*shat*shat*shat*std::log(shat)*std::log(shat)+18.*PI*I*shat*(1.-2.*shat)*std::log(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)+PI2*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
}

complex_t BKstarllDecay::B_Seidel(double s, double m_b, double mu_b) {
    double shat = s/m_b/m_b;
	double z = 4./shat;

	complex_t x1=0.5+0.5*I*std::sqrt(z-1.);
	complex_t x2=0.5-0.5*I*std::sqrt(z-1.);
	complex_t x3=0.5+0.5*I/std::sqrt(z-1.);
	complex_t x4=0.5-0.5*I/std::sqrt(z-1.);

	if (fabs(shat - 1.0) < 1e-2) {
		complex_t mu_terms =
		8./243./shat * (4.-34.*shat-17.*PI*I*shat) * log(m_b*m_b/mu_b/mu_b)
		+ 17.*shat * 8./243./shat * log(shat) * log(m_b*m_b/mu_b/mu_b)
		+ (2.+shat) * std::sqrt(z-1.)/729./shat * (-48.) * (PI/2. - std::atan(std::sqrt(z-1.))) * log(m_b*m_b/mu_b/mu_b);

		complex_t nonmu_terms = -1.2534705628994441 + 3.1545210184193809 * I +
		(-1.1399966466176837 - 1.3704066719362884 * I) * (shat - 1.0) +
		(0.77575942579740349 + 0.59987612809286587 * I) * pow(shat - 1.0, 2);

		complex_t result = mu_terms + nonmu_terms;
		return result;
	}

	/* In the A expression, for (CLi2(shat)+log(shat)*std::log(1.-shat)), the real part is calculated correctly
	 * but the imaginary part above the branch cut is not (in principle it should get cancelled).
	 * I fixed the CLi2 to use the same branch cut as std::log, so it can be used directly */
	complex_t Li2log_term = CLi2(shat)+std::log(shat)*std::log(1.-shat);

	return 8./243./shat*((4.-34.*shat-17.*PI*I*shat)*log(m_b*m_b/mu_b/mu_b)+8.*shat*pow(log(m_b*m_b/mu_b/mu_b),2.)+17.*shat*std::log(shat)*log(m_b*m_b/mu_b/mu_b))
	+(2.+shat)*std::sqrt(z-1.)/729./shat*(-48.*log(m_b*m_b/mu_b/mu_b)*(PI/2.-std::atan(std::sqrt(z-1.)))-18.*PI*std::log(z-1.)+3.*I*std::log(z-1.)*std::log(z-1.)
	-24.*I*CLi2(-x2/x1)-5.*PI2*I+6.*I*(-9.*std::log(x1)*std::log(x1)+std::log(x2)*std::log(x2)-2.*std::log(x4)*std::log(x4)+6.*std::log(x1)*std::log(x2)-4.*std::log(x1)*std::log(x3)+8.*std::log(x1)*std::log(x4))
	-12.*PI*(2.*std::log(x1)+std::log(x3)+std::log(x4)))
	-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*Li2log_term
	+3.*(2.+shat)*(3.-shat)*std::log(x2/x1)*std::log(x2/x1)+12.*PI*(-6.-shat+shat*shat)*(PI/2.-std::atan(std::sqrt(z-1.))))
	+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*std::log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*PI*I*shat*(82.-173.*shat+73.*shat*shat))
	-4./243./shat/pow(1.-shat,3.)*(-2.*std::sqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(PI/2.-std::atan(std::sqrt(z-1.)))-9.*shat*shat*shat*std::log(shat)*std::log(shat)+2.*PI*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*std::log(shat))
	+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)-PI2*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
}

complex_t BKstarllDecay::C_Seidel(double s, double mu_b) {
    return -16./81.*std::log(s/mu_b/mu_b)+428./243.-64./27.*ZETA3+16./81.*PI*I;
}

double BKstarllDecay::F_a(FF a, double q2) {
    if (!FF_loaded) {
        LOG_ERROR("Logic Error", "Form Factor parameters were used before initialization.");
    }

    auto ai = cache.alpha_ai.at(a);
    double P = pole(q2, cache.m_R.at(a));
    double Z = z(q2) - cache.z0;

    return P * (ai[0] + Z * (ai[1] + Z * ai[2]));
}

double BKstarllDecay::E_K(double q2, double m_B, double m_K) {
    return (std::pow(m_B, 2) + std::pow(m_K, 2) - q2) / (2 * m_B);
}

double BKstarllDecay::A_2(double q2, double m_B, double m_K) {
    double A_1 = F_a(FF::A1, q2);
    double A_12 = F_a(FF::A12, q2);
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    return (cache.tp * (mB2 - mK2 - q2) * A_1 - 16 * m_B * mK2 * (m_B + m_K) * A_12) / ((cache.tp - q2) * (cache.tm - q2));
}

double BKstarllDecay::T_3(double q2, double m_B, double m_K) {
    double T_2 = F_a(FF::T2, q2);
    double T_23 = F_a(FF::T23, q2);
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    return ((mB2 - mK2) * (mB2 + 3 * mK2 - q2) * T_2 - 8 * m_B * mK2 * (m_B - m_K) * T_23) / ((cache.tp - q2) * (cache.tm - q2));
}

double BKstarllDecay::xi_perp(double q2, double m_B, double m_K) {
    return m_B * F_a(FF::V, q2) / (m_B + m_K);
}

double BKstarllDecay::xi_par(double q2, double m_B, double m_K) {
    return (m_B + m_K) * F_a(FF::A1, q2) / (2 * E_K(q2, m_B, m_K)) - (m_B - m_K) * A_2(q2, m_B, m_K) / m_B;
}

complex_t BKstarllDecay::T_par_m_0(double m_B) {
    int delta_qu = cfg.charge == BKstarllConfig::B_Charge::B_PLUS;
    return 4. * m_B / cache.m_b_PS * (delta_qu * 3. * cache.lambda_hat_u * cache.C[WCoef::C2] - cache.C_bar[WCoef::C3] - 3. * cache.C_bar[WCoef::C4]);
}

complex_t BKstarllDecay::T_par_p_p_f(double u, double q2, double m_B, double m_K) {
    return 2. * T_perp_p_p_f(u, q2, m_B, m_K);
}

complex_t BKstarllDecay::T_par_p_m_f(double u, double q2, double m_B, double m_K) {
    return 2. * T_perp_p_m_f(u, q2, m_B, m_K);
}

complex_t BKstarllDecay::T_perp_p_p_f(double u, double q2, double m_B, double m_K) {
    return 2. * m_B / (1. - u) / E_K(q2, m_B, m_K) * (cache.C[WCoef::C7] + cache.C[WCoef::CP7]);
}

complex_t BKstarllDecay::T_perp_p_m_f(double u, double q2, double m_B, double m_K) {
    return 2. * m_B / (1. - u) / E_K(q2, m_B, m_K) * (cache.C[WCoef::C7] - cache.C[WCoef::CP7]);
}

complex_t BKstarllDecay::T_perp_p_nf(double u, double q2, double m_B, double m_K) {
    double E = E_K(q2, m_B, m_K);
    complex_t t_perp_mc = t_perp(u, cache.m_c_pole, q2, E, m_B);
    complex_t t_perp_mb = t_perp(u, cache.m_b_PS, q2, E, m_B);
    complex_t t_perp_0 = t_perp(u, 0, q2, E, m_B);
    return -4 * cache.e_d * cache.C[WCoef::C8] / (u + (1 - u) * q2 / (m_B * m_B))
            + m_B / (2 * cache.m_b_PS) * (
                cache.e_u * (
                    t_perp_mc * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] + cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.))
                    - t_perp_0 * cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)
                )
                + cache.e_d * (
                    t_perp_mb * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] - 4 * cache.m_b_PS / m_B * cache.C_bar[WCoef::C5])
                    + t_perp_0 * cache.C_bar[WCoef::C3]
                )
            );
}

complex_t BKstarllDecay::T_par_p_nf(double u, double q2, double m_B, double m_K) {
    double E = E_K(q2, m_B, m_K);
    complex_t t_par_mc = t_par(u, cache.m_c_pole, q2, E, m_B);
    complex_t t_par_mb = t_par(u, cache.m_b_PS, q2, E, m_B);
    complex_t t_par_0 = t_par(u, 0, q2, E, m_B);
    return m_B / cache.m_b_PS * (
        cache.e_u * (
            t_par_mc * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] + cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.))
            - t_par_0 * cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)
        )
        + cache.e_d * (
            t_par_mb * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6])
            + t_par_0 * cache.C_bar[WCoef::C3]
        )
    );
}

complex_t BKstarllDecay::T_par_m_nf(double u, double q2, double m_B, double m_K) {
    double E = E_K(q2, m_B, m_K);
    double v = m_B * m_B * (1 - u) + q2 * u;
    return 8. * m_B * m_B * cache.C[WCoef::C8] / v + 6 * m_B / cache.m_b_PS * (
        h(v, cache.m_c_pole, w_config.hadronic_scale) * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6] + cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
        + h(v, cache.m_b_pole, w_config.hadronic_scale) * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6]) 
        + h(v, 0, w_config.hadronic_scale) * (cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4] + 3. * cache.C_bar[WCoef::C6] - cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
        - 8. / 27. * (cache.C_bar[WCoef::C3] - cache.C_bar[WCoef::C5] - 15. * cache.C_bar[WCoef::C6])
    );
}

complex_t BKstarllDecay::inv_lambda_B_m(double q2, double m_B) {
    double x = q2 / (m_B * cache.omega_0);
    return std::exp(-x) / cache.omega_0 * (I * PI - Ei(x));
}

complex_t BKstarllDecay::I_perp_p(double q2, double m_B, double m_K) {
    auto f = [q2, m_B, m_K, this] (double u) {
        return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * (T_perp_p_p_f(u, q2, m_B, m_K) + T_perp_p_nf(u, q2, m_B, m_K));
    };
    
    return cache.alpha_s_mu_f / (4 * PI) * cache.C_F / cache.lambda_B_p * c_integrate(f, 0, 1, 1e-3);
}

complex_t BKstarllDecay::I_perp_m(double q2, double m_B, double m_K) {
    auto f = [q2, m_B, m_K, this] (double u) {
        return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * (T_perp_p_m_f(u, q2, m_B, m_K) + T_perp_p_nf(u, q2, m_B, m_K));
    };
    
    return cache.alpha_s_mu_f / (4 * PI) * cache.C_F / cache.lambda_B_p * c_integrate(f, 0, 1, 1e-3);
}

complex_t BKstarllDecay::I_par_p(double q2, double m_B, double m_K) {
    auto f = [q2, m_B, m_K, this] (double u) {
        double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        complex_t i1 = phi * (T_par_p_p_f(u, q2, m_B, m_K) + T_par_p_nf(u, q2, m_B, m_K));
        complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K));
        return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2, m_B) * i2;
    };
    
    return c_integrate(f, 0, 1, 1e-3);
}

complex_t BKstarllDecay::I_par_m(double q2, double m_B, double m_K) {
    auto f = [q2, m_B, m_K, this] (double u) {
        double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        complex_t i1 = phi * (T_par_p_m_f(u, q2, m_B, m_K) + T_par_p_nf(u, q2, m_B, m_K));
        complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K));
        return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2, m_B) * i2;
    };
    
    return c_integrate(f, 0, 1, 1e-3);
}


void BKstarllDecay::build_op_tree() {

    // SM Parameters
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));

    // Wilson node
    auto wilson = this->get_wilson_node();

    // Formfactor node
    auto ff = this->load_FF_params();
    
    // Flavor Parameters
    auto m_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));
    auto m_Ks = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 323));

    auto tp = std::make_shared<OperatorNode>("t_+", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.tp = pow(values[0] + values[1], 2); return cache.tp; });
    tp->addChildren({m_B, m_Ks});
    auto tm = std::make_shared<OperatorNode>("t_-", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.tm = pow(values[0] - values[1], 2); return cache.tm; });
    tm->addChildren({m_B, m_Ks});
    auto t0 = std::make_shared<OperatorNode>("t_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        cache.t0 = this->cfg.ff_src == BKstarllConfig::FF_Src::HLMW ? 12. : cache.tp * (1. - std::sqrt(1 - cache.tm / cache.tp));
        return cache.t0; 
    });
    t0->addChildren({tp, tm});
    auto z0 = std::make_shared<OperatorNode>("z_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return z(0); });
    z0->addChildren({tp, t0});

    auto test_ff = std::make_shared<OperatorNode>("test_ff", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {

        std::ofstream fs;
        fs.open("B_Ksll_FF.csv");
        fs << "q2,A0,A1,A12,V,T1,T2,T23\n";

        auto write_line = [&] (double q2) {
            fs << q2 << "," 
            << F_a(FF::A0, q2) << ","
            << F_a(FF::A1, q2) << ","
            << F_a(FF::A12, q2) << ","
            << F_a(FF::V, q2) << ","
            << F_a(FF::T1, q2) << ","
            << F_a(FF::T2, q2) << ","
            << F_a(FF::T23, q2)
            << "\n";
        };

        double q2_min = 4 * std::pow(0.1257, 2);
        double q2_max = 19.2542;
        size_t n = 200;
        double dq2 = (q2_max - q2_min) / n;
        double q2 = q2_min;
        for (size_t i = 0; i <= n; i++) {
            write_line(q2);
            q2 += dq2;
        }

        return 0; 
    });
    test_ff->addChildren({tp, t0, z0, ff});
    roots.emplace(Observables::TEST_B__KS_L_L, test_ff);
}