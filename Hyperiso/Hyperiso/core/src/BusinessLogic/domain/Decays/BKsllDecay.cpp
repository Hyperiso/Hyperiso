#include "BKsllDecay.h"

using Charge = BKstarllConfig::B_Charge;
using FF_Src = BKstarllConfig::FF_Src;

std::shared_ptr<OperatorNode> BKstarllDecay::load_FF_params() {
    auto ff_params = std::make_shared<OperatorNode>("FF_params", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 0; });
    int ff_id = (int)(this->cfg.ff_src) + 1;
    int sse_order = this->cfg.ff_src == FF_Src::HLMW ? 1 : 2;
    
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
            - 0.5 * h(q2, 0., w_config.hadronic_scale) * (cache.C[WCoef::C3] + 4./3. * cache.C[WCoef::C4] + 16. * cache.C[WCoef::C5] + 64./3. * cache.C[WCoef::C6])
            + 4./3. * cache.C[WCoef::C3] + 64./9. * cache.C[WCoef::C5] + 64./27. * cache.C[WCoef::C6];
}

complex_t BKstarllDecay::Y_u(double q2) {
    return (h(q2, cache.m_c_pole, w_config.hadronic_scale) - h(q2, 0, w_config.hadronic_scale)) * (4. / 3. * cache.C[WCoef::C1] + cache.C[WCoef::C2]);
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

void BKstarllDecay::fill_wilson_cache() {
    auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    auto bq_wilsons = w_proxy->getAFR(WGroup::BScalar, this->w_config.order);
    WCoef bp_cached[5] {WCoef::CP7, WCoef::CP9, WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2};

    for (auto p : b_wilsons) cache.C.emplace(p); 
    for (auto p : bq_wilsons) cache.C.emplace(p);
    for (auto id : bp_cached) cache.C.emplace(std::pair{id, bp_wilsons.at(id)});
}

void BKstarllDecay::fill_wilson_bar_cache() {
    auto b_ids = WCoefMapper::B_group();

    for (size_t i = 0; i < 6; i++) {
        cache.C_bar[b_ids[i]] = 0;
        for (size_t j = 0; j < 6; j++) {
            cache.C_bar[b_ids[i]] += P_bar[i][j] * cache.C[b_ids[j]];
        }
    }
}

std::shared_ptr<ParameterNode> BKstarllDecay::get_lepton_mass() {
    switch (cfg.gen) {
        case BKstarllConfig::Lepton::E:
            return std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 11));
        case BKstarllConfig::Lepton::MU:
            return std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 13));
        case BKstarllConfig::Lepton::TAU:
            return std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 15));
        default:
            return nullptr;
    }
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

double BKstarllDecay::F_perp(double s) {
    double d = s - 1;
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double s2 = s * s;
    double ls = std::log(s);
    double f0 = (s + 1) / d2 - 2 * s * ls / d3;
    double f1 = -(s2 + 10 * s + 1) / d3 + 6 * s * (s + 1) * ls / d4;
    double f2 = (s + 1) + (s2 + 28 * s + 1) / d4 - 12 * s * (s2 + 3 * s + 1) * ls / d5;
    return f0 + cache.a_1_perp * f1 + cache.a_2_perp * f2;
}

double BKstarllDecay::X_perp(double s) {
    double d = s - 1;
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double s2 = s * s;
    double ls = std::log(s);
    double f0 = (s - 3) / d2 + 2 * ls / d3;
    double f1 = -(s2 - 8 * s - 17) / d3 + 6 * (3 * s + 1) * ls / d4;
    double f2 = -(s * s2 - 15 * s2 - 123 * s - 43) / d4 - 12 * (6 * s2 + 8 * s + 1) * ls / d5;
    return f0 + cache.a_1_perp * f1 + cache.a_2_perp * f2;
}

double BKstarllDecay::gv_dga_4(double u, double z3a, double z3v, double w10a, double dtp, double dtm) {
    double a1 = -60 * z3a * (w10a + 4) + 1680 * z3v;
    double a2 = 30 * z3a * (15 * w10a + 32) - 12600 * z3v + 36 * cache.a_1_par - 72 * cache.a_2_par - 12;
    double a3 = -100 * z3a * (9 * w10a + 8) + 25200 * z3v - 48 * cache.a_1_par + 240 * cache.a_2_par;
    double a4 = 525 * z3a * w10a - 14700 * z3v - 180 * cache.a_2_par;
    return -u * (a1 + u * (a2 + u * (a3 + u * a4))) / 4 + dtp * (9 * u - 1.5) + dtm * 6 * u + 3 * (dtp + dtm) * log(1 - u);
}

complex_t BKstarllDecay::F_V(double v) {
    return .75 * (
        h(v, cache.m_c_pole, w_config.hadronic_scale) * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6] + cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
      + h(v, cache.m_b_pole, w_config.hadronic_scale) * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6]) 
      + h(v, 0, w_config.hadronic_scale) * (cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4] + 3. * cache.C_bar[WCoef::C6] - cache.lambda_hat_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
      - 8. / 27. * (cache.C_bar[WCoef::C3] - cache.C_bar[WCoef::C5] - 15. * cache.C_bar[WCoef::C6])
    );
}

double BKstarllDecay::L(double q2) {
    double mb2 = std::pow(cache.m_b_PS, 2);
    return (q2 - mb2) * std::log(1 - q2 / mb2) / q2;
}

complex_t BKstarllDecay::C_perp_0(double q2, double m_B, double sign) {
    return (cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7]) + q2 * (Y(q2) + cache.lambda_hat_u * Y_u(q2)) / (2 * cache.m_b_PS * m_B);
}

complex_t BKstarllDecay::C_par_0(double q2, double m_B, double sign) {
    return -(cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7]) - m_B * (Y(q2) + cache.lambda_hat_u * Y_u(q2)) / (2 * cache.m_b_PS);
}

complex_t BKstarllDecay::C_perp_f(double q2, double sign) {
    return (cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7]) * (2 * std::log(cache.m_b_PS / w_config.hadronic_scale) - L(q2) + cache.Delta_M);
}

complex_t BKstarllDecay::C_par_f(double q2, double sign) {
    return -(cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7]) * (2 * std::log(cache.m_b_PS / w_config.hadronic_scale) + 2 * L(q2) + cache.Delta_M);
}

complex_t BKstarllDecay::C_perp_nf(double q2, double m_B) {
    double s_hat = q2 / (cache.m_b_PS * cache.m_b_PS);
    complex_t F_27 = f_27(s_hat, cache.L_b, cache.z_c) * (1. + cache.lambda_hat_u) + F_27_u(s_hat, cache.L_b) * cache.lambda_hat_u;
    complex_t F_19 = f_19_PS(s_hat, cache.L_b, cache.z_c) * (1. + cache.lambda_hat_u) + F_19_u(s_hat, cache.L_b) * cache.lambda_hat_u;
    complex_t F_29 = f_29_PS(s_hat, cache.L_b, cache.z_c) * (1. + cache.lambda_hat_u) + F_29_u(s_hat, cache.L_b) * cache.lambda_hat_u;
    return -(
        cache.C_bar[WCoef::C2] * F_27 
      + cache.C[WCoef::C8] * f_87(s_hat, cache.L_b)
      + q2 / (2 * cache.m_b_PS * m_B) * (
            (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C1] / 3.) * F_29
          + 2. * cache.C_bar[WCoef::C1] * F_19
          + cache.C[WCoef::C8] * f_89(s_hat)
        )
    ) / cache.C_F;
}

complex_t BKstarllDecay::C_par_nf(double q2, double m_B) {
    double s_hat = q2 / (m_B * m_B);
    complex_t F_27 = f_27(s_hat, cache.L_b, cache.z_c) * (1. + cache.lambda_hat_u) + F_27_u(s_hat, cache.L_b) * cache.lambda_hat_u;
    complex_t F_19 = f_19_PS(s_hat, cache.L_b, cache.z_c) * (1. + cache.lambda_hat_u) + F_19_u(s_hat, cache.L_b) * cache.lambda_hat_u;
    complex_t F_29 = f_29_PS(s_hat, cache.L_b, cache.z_c) * (1. + cache.lambda_hat_u) + F_29_u(s_hat, cache.L_b) * cache.lambda_hat_u;
    return (
        cache.C_bar[WCoef::C2] * F_27
      + cache.C[WCoef::C8] * f_87(s_hat, cache.L_b)
      + m_B / (2 * cache.m_b_PS) * (
            (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C1] / 3.) * F_29
           + 2. * cache.C_bar[WCoef::C1] * F_19
           + cache.C[WCoef::C8] * f_89(s_hat))
    ) / cache.C_F;
}

complex_t BKstarllDecay::T_par_m_0(double m_B) {
    int delta_qu = cfg.charge == Charge::B_PLUS;
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
    double v = m_B * m_B * (1 - u) + q2 * u;
    return 8. * m_B * m_B * cache.C[WCoef::C8] / v + 8 * m_B / cache.m_b_PS * F_V(v);
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
    
    return c_integrate(f, 0, 1 - 1e-6, 1e-3);
}

complex_t BKstarllDecay::I_par_m(double q2, double m_B, double m_K) {
    auto f = [q2, m_B, m_K, this] (double u) {
        double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        complex_t i1 = phi * (T_par_p_m_f(u, q2, m_B, m_K) + T_par_p_nf(u, q2, m_B, m_K));
        complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K));
        return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2, m_B) * i2;
    };
    
    return c_integrate(f, 0, 1 - 1e-6, 1e-3);
}

complex_t BKstarllDecay::I_HSA_1(double q2, double m_B) {
    auto f = [q2, m_B, this] (double u) {
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        double v = m_B * m_B * (1 - u) + u * q2;
        return phi * m_B * m_B / v * F_V(v);
    };
    
    return c_integrate(f, 0, 1, 1e-3);
}

complex_t BKstarllDecay::I_HSA_2(double q2, double m_B, double z3a, double z3v, double w10a, double dtp, double dtm) {
    auto f = [q2, m_B, z3a, z3v, w10a, dtp, dtm, this] (double u) {
        double int_phi_par = gv_dga_4(u, z3a, z3v, w10a, dtp, dtm);
        double v = m_B * m_B * (1 - u) + u * q2;
        return int_phi_par * F_V(v);
    };
    
    return c_integrate(f, 0, 1, 1e-3);
}

complex_t BKstarllDecay::delta_T_perp_WA(double q2, double m_B, double m_K, double f_B, double f_K_par) {
    double pref = cache.e_q * 2 * PI2 * f_B / (cache.m_b_PS * m_B);
    complex_t W_perp = cache.C[WCoef::C3] + 4. / .3 * (cache.C[WCoef::C4] + 3. * cache.C[WCoef::C5] + 4. * cache.C[WCoef::C6]);
    complex_t W_par = cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4];
    if (cfg.charge == Charge::B_PLUS) 
        W_par += -3. * cache.C[WCoef::C2];
    double s_hat = q2 / (m_B * m_B);
    return pref * (
        -2 * cache.f_K_perp * W_perp * F_perp(s_hat)
      + f_K_par * m_K * W_par / (3 * (1 - s_hat) * cache.lambda_B_p)
    );
}

complex_t BKstarllDecay::delta_T_perp_HSA(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm) {
    double pref = cache.e_q * cache.alpha_s_mu_b * cache.C_F * PI * f_B / (cache.Nc * cache.m_b_PS * m_B);
    double s_hat = q2 / (m_B * m_B);
    return pref * (
        3. * cache.C[WCoef::C8] * cache.m_b_PS / m_B * cache.f_K_perp * X_perp(s_hat)
      + 2 * cache.f_K_perp * I_HSA_1(q2, m_B)
      - m_K * f_K_par / ((1 * s_hat) * cache.lambda_B_p) * I_HSA_2(q2, m_B, z3a, z3v, w10a, dtp, dtm)
    );
}

complex_t BKstarllDecay::T_perp_p(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm) {
    complex_t C_perp_p = C_perp_0(q2, m_B, 1) + cache.alpha_s_mu_b / (4 * PI) * (C_perp_f(q2, 1) + C_perp_nf(q2, m_B));
    return xi_perp(q2, m_B, m_K) * C_perp_p + cache.pref_T_perp * I_perp_p(q2, m_B, m_K) + delta_T_perp_WA(q2, m_B, m_K, f_B, f_K_par) + delta_T_perp_HSA(q2, m_B, m_K, f_B, f_K_par, z3a, z3v, w10a, dtp, dtm);
}

complex_t BKstarllDecay::T_perp_m(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm) {
    complex_t C_perp_m = C_perp_0(q2, m_B, -1) + cache.alpha_s_mu_b / (4 * PI) * (C_perp_f(q2, -1) + C_perp_nf(q2, m_B));
    return xi_perp(q2, m_B, m_K) * C_perp_m + cache.pref_T_perp * I_perp_m(q2, m_B, m_K) + delta_T_perp_WA(q2, m_B, m_K, f_B, f_K_par) + delta_T_perp_HSA(q2, m_B, m_K, f_B, f_K_par, z3a, z3v, w10a, dtp, dtm);
}

complex_t BKstarllDecay::T_par_p(double q2, double m_B, double m_K) {
    complex_t C_par_p = C_par_0(q2, m_B, 1) + cache.alpha_s_mu_b / (4 * PI) * (C_par_f(q2, 1) + C_par_nf(q2, m_B));
    return xi_par(q2, m_B, m_K) * C_par_p + cache.pref_T_par / E_K(q2, m_B, m_K) * I_par_p(q2, m_B, m_K);
}

complex_t BKstarllDecay::T_par_m(double q2, double m_B, double m_K) {
    complex_t C_par_m = C_par_0(q2, m_B, -1) + cache.alpha_s_mu_b / (4 * PI) * (C_par_f(q2, -1) + C_par_nf(q2, m_B));
    return xi_par(q2, m_B, m_K) * C_par_m + cache.pref_T_par / E_K(q2, m_B, m_K) * I_par_m(q2, m_B, m_K);
}

complex_t BKstarllDecay::Delta_par(double q2, double m_B, double m_K, double f_B, double f_K_par) {
    return 1 + cache.alpha_s_mu_b * cache.C_F / (2 * PI) * (
        L(q2) - 1 - 3 * PI2 * q2 * f_B * f_K_par * m_K / (cache.Nc * m_B * cache.lambda_B_p * xi_par(q2, m_B, m_K) * std::pow(E_K(q2, m_B, m_K), 3)) * F_perp(0.0)
    );
}

complex_t BKstarllDecay::T_perp_p_cached(double q2) {
    return lerp(q2, T_perp_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_perp_m_cached(double q2) {
    return lerp(q2, T_perp_m_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_par_p_cached(double q2) {
    return lerp(q2, T_par_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_par_m_cached(double q2) {
    return lerp(q2, T_par_m_lookup, cache.q2_min, cache.q2_high);
}

double BKstarllDecay::beta_l(double q2, double m_l) {
    return std::sqrt(1 - 4 * m_l * m_l / q2);
}

double BKstarllDecay::lambda(double q2, double m_B, double m_K) {
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    return mB2 * mB2 + mK2 * mK2 + q2 * q2 - 2 * (mB2 * mK2 + (mB2 + mK2) * q2);
}

complex_t BKstarllDecay::N(double q2, double m_B, double m_K, double m_l) {
    return std::sqrt(q2 * beta_l(q2, m_l) * std::sqrt(lambda(q2, m_B, m_K)));
}

complex_t BKstarllDecay::A_perp(double q2, double m_B, double m_K, double m_l, double sign) {
    return N(q2, m_B, m_K, m_l) * std::sqrt(2 * lambda(q2, m_B, m_K)) * (
        (cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10])) * F_a(FF::V, q2) / (m_B + m_K)
      + 2 * cache.m_b_PS * T_perp_p_cached(q2) / q2
    );
}

complex_t BKstarllDecay::A_par(double q2, double m_B, double m_K, double m_l, double sign) {
    return -N(q2, m_B, m_K, m_l) * std::sqrt(2.) * (m_B * m_B - m_K * m_K) * (
        (cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10])) * F_a(FF::A1, q2) / (m_B - m_K)
      + 4. * cache.m_b_PS * E_K(q2, m_B, m_K) * T_perp_m_cached(q2) / (m_B * q2)
    );
}

complex_t BKstarllDecay::A_0(double q2, double m_B, double m_K, double m_l, double sign) {
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    return -N(q2, m_B, m_K, m_l) / (2. * m_K * std::sqrt(q2)) * (
        (cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10])) * (
            (mB2 - mK2 - q2) * (m_B + m_K) * F_a(FF::A1, q2)
          - lambda(q2, m_B, m_K) * A_2(q2, m_B, m_K) / (m_B + m_K)
        )
      + 2. * cache.m_b_PS * (
          2. * E_K(q2, m_B, m_K) * (mB2 + 3. * mK2 - q2) * T_perp_m_cached(q2) 
        - lambda(q2, m_B, m_K) * (T_perp_m_cached(q2) + T_par_m_cached(q2)) / (mB2 - mK2)
      )
    );
}

complex_t BKstarllDecay::A_t(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    return N(q2, m_B, m_K, m_l) * std::sqrt(lambda(q2, m_B, m_K) / q2) * (
        (cache.C[WCoef::C10] - cache.C[WCoef::CP10] + q2 / (m_l * (cache.m_b_mu_b + m_s)) * (cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2]))
    ) * E_K(q2, m_B, m_K) * xi_par(q2, m_B, m_K) / (m_K * Delta_par(q2, m_B, m_K, f_B, f_K_par));
}

complex_t BKstarllDecay::A_S(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    return -2. * N(q2, m_B, m_K, m_l) / (cache.m_b_mu_b + m_s) * std::sqrt(lambda(q2, m_B, m_K)) 
            * (cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1]) 
            * E_K(q2, m_B, m_K) * xi_par(q2, m_B, m_K) / (m_K * Delta_par(q2, m_B, m_K, f_B, f_K_par));
}

double BKstarllDecay::J1s(double q2, double m_B, double m_K, double m_l) {
    return (2. + std::pow(beta_l(q2, m_l), 2)) / 4. * (
        std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1)), 2) 
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1)), 2)
      + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1)), 2)
      + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1)), 2)
    ) + std::pow(2 * m_l, 2) / q2 * std::real(
        A_perp(q2, m_B, m_K, m_l, -1) * std::conj(A_perp(q2, m_B, m_K, m_l, 1))
      + A_par(q2, m_B, m_K, m_l, -1) * std::conj(A_par(q2, m_B, m_K, m_l, 1))
    );
}

double BKstarllDecay::J1c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    return std::pow(std::abs(A_0(q2, m_B, m_K, m_l, -1)), 2) + std::pow(std::abs(A_0(q2, m_B, m_K, m_l, 1)), 2) 
         + std::pow(2 * m_l, 2) / q2 * (
              std::pow(std::abs(A_t(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)), 2)
            + 2. * std::real(A_0(q2, m_B, m_K, m_l, -1) * std::conj(A_0(q2, m_B, m_K, m_l, 1)))
           )
         + std::pow(beta_l(q2, m_l) * std::abs(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)), 2);
}

double BKstarllDecay::J2s(double q2, double m_B, double m_K, double m_l) {
    return std::pow(beta_l(q2, m_l), 2) / 4. * (
        std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1)), 2) 
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1)), 2)
      + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1)), 2)
      + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1)), 2)
    );
}

double BKstarllDecay::J2c(double q2, double m_B, double m_K, double m_l) {
    return -std::pow(beta_l(q2, m_l), 2) * (
        std::pow(std::abs(A_0(q2, m_B, m_K, m_l, -1)), 2) 
      + std::pow(std::abs(A_0(q2, m_B, m_K, m_l, 1)), 2)
    );
}

double BKstarllDecay::J3(double q2, double m_B, double m_K, double m_l) {
    return std::pow(beta_l(q2, m_l), 2) / 2. * (
        std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1)), 2) 
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1)), 2)
      - std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1)), 2)
      - std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1)), 2)
    );
}

double BKstarllDecay::J4(double q2, double m_B, double m_K, double m_l) {
    return std::pow(beta_l(q2, m_l), 2) / std::sqrt(2.) * (
        std::real(A_0(q2, m_B, m_K, m_l, -1) * std::conj(A_par(q2, m_B, m_K, m_l, -1))) 
      + std::real(A_0(q2, m_B, m_K, m_l, 1) * std::conj(A_par(q2, m_B, m_K, m_l, 1)))
    );
}

double BKstarllDecay::J5(double q2, double m_B,double m_K,double m_l,double f_B,double f_K_par,double m_s) {
    return beta_l(q2, m_l) * std::sqrt(2.) * (
        std::real(A_0(q2, m_B, m_K, m_l, -1) * std::conj(A_perp(q2, m_B, m_K, m_l, -1))) 
      - std::real(A_0(q2, m_B, m_K, m_l, 1) * std::conj(A_perp(q2, m_B, m_K, m_l, 1)))
      - m_l / std::sqrt(q2) * std::real((A_par(q2, m_B, m_K, m_l, -1) + A_par(q2, m_B, m_K, m_l, 1)) * std::conj(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)))
    );
}

double BKstarllDecay::J6s(double q2, double m_B, double m_K, double m_l) {
    return 2. * beta_l(q2, m_l) * (
        std::real(A_par(q2, m_B, m_K, m_l, -1) * std::conj(A_perp(q2, m_B, m_K, m_l, -1))) 
      - std::real(A_par(q2, m_B, m_K, m_l, 1) * std::conj(A_perp(q2, m_B, m_K, m_l, 1))) 
    );
}

double BKstarllDecay::J6c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    return 4. * beta_l(q2, m_l) * m_l / std::sqrt(q2) * (std::real((A_0(q2, m_B, m_K, m_l, -1) + A_0(q2, m_B, m_K, m_l, 1)) * std::conj(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s))));
}

double BKstarllDecay::J7(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    return beta_l(q2, m_l) * std::sqrt(2.) * (
        std::imag(A_0(q2, m_B, m_K, m_l, -1) * std::conj(A_par(q2, m_B, m_K, m_l, -1))) 
      - std::imag(A_0(q2, m_B, m_K, m_l, 1) * std::conj(A_par(q2, m_B, m_K, m_l, 1)))
      + m_l / std::sqrt(q2) * std::imag((A_perp(q2, m_B, m_K, m_l, -1) + A_perp(q2, m_B, m_K, m_l, 1)) * std::conj(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)))
    );
}

double BKstarllDecay::J8(double q2, double m_B, double m_K, double m_l) {
    return std::pow(beta_l(q2, m_l), 2) / std::sqrt(2.) * (
        std::imag(A_0(q2, m_B, m_K, m_l, -1) * std::conj(A_perp(q2, m_B, m_K, m_l, -1))) 
      + std::imag(A_0(q2, m_B, m_K, m_l, 1) * std::conj(A_perp(q2, m_B, m_K, m_l, 1)))
    );
}

double BKstarllDecay::J9(double q2, double m_B, double m_K, double m_l) {
    return std::pow(beta_l(q2, m_l), 2) * (
        std::imag(A_par(q2, m_B, m_K, m_l, -1) * std::conj(A_perp(q2, m_B, m_K, m_l, -1))) 
      + std::imag(A_par(q2, m_B, m_K, m_l, 1) * std::conj(A_perp(q2, m_B, m_K, m_l, 1))) 
    );
}


void BKstarllDecay::build_op_tree() {

    // SM Parameters
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));
    auto inv_alpha_em = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 1));
    auto m_b_pole = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "QCD", {5, 2}));
    auto m_c_pole = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "QCD", 4));
    auto V_ub = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {0, 2}));
    auto V_us = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {0, 1}));
    auto V_tb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {2, 2}));
    auto V_ts = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", {2, 1}));
    auto m_l = get_lepton_mass();
    auto m_s = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 3));

    // Wilson node
    auto wilson = this->get_wilson_node();
    auto wilson_cache = std::make_shared<OperatorNode>("wilson_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { fill_wilson_cache(); return 0; });
    wilson_cache->addChild(wilson);
    auto wilson_bar_cache = std::make_shared<OperatorNode>("wilson_bar_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { fill_wilson_bar_cache(); return 0; });
    wilson_bar_cache->addChild(wilson_cache);

    // Formfactor node
    auto ff = this->load_FF_params();
    
    // Flavor Parameters
    // TODO : Manage the CHARGE flag (B0 - K*0 // B+ - K*+)
    auto m_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));
    auto m_Ks = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 323));
    auto f_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", {511, 1}));
    auto f_K_perp_1GeV = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", {323, 1}));
    auto f_K_par = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", {323, 2}));

    // Decay-specific parameters
    auto a_1_perp_1GeV = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 1));
    auto a_2_perp_1GeV = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 2));
    auto a_1_par_1GeV = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 3));
    auto a_2_par_1GeV = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 4));
    auto zeta_3_A = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 5));
    auto zeta_3_V = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 6));
    auto omega_10_A = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 7));
    auto delta_tilde_p = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 8));
    auto delta_tilde_m = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 9));
    auto lambda_B_1GeV = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 10));
    auto lambda_h = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ks", 12));
    auto q2_low = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ksll", {10, 1}));
    auto q2_high = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Ksll", {10, 2}));

    auto n_q2_low = std::make_shared<OperatorNode>("q2_low", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.q2_low = values[0]; return cache.q2_low; });
    n_q2_low->addChild(q2_low);
    auto n_q2_high = std::make_shared<OperatorNode>("q2_high", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.q2_high = values[0]; return cache.q2_high; });
    n_q2_high->addChild(q2_high);
    auto q2_min = std::make_shared<OperatorNode>("q2_min", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.q2_min = pow(2. * values[0], 2); return cache.q2_min; });
    q2_min->addChildren({m_l});
    auto q2_max = std::make_shared<OperatorNode>("q2_max", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.q2_max = pow(values[0] - values[1], 2); return cache.q2_max; });
    q2_max->addChildren({m_B, m_Ks});
    auto eq = std::make_shared<OperatorNode>("e_q", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.e_q = cfg.charge == Charge::B_0 ? cache.e_d : cache.e_u; return cache.e_q; });
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
    auto C_F = std::make_shared<OperatorNode>("C_F", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.C_F = ObsQCDProxy().get_constants()->C_F; return cache.C_F; });
    auto N_c = std::make_shared<OperatorNode>("N_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.Nc = ObsQCDProxy().get_constants()->Nc; return cache.Nc; });
    auto beta_0 = std::make_shared<OperatorNode>("beta_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.beta_0 = ObsQCDProxy().get_constants()->beta[0][QCDHelper::get_nf(w_config.hadronic_scale)]; return cache.beta_0; });
    auto mu_f = std::make_shared<OperatorNode>("mu_f", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.mu_f = sqrt(values[0] * w_config.hadronic_scale); return cache.mu_f; });
    mu_f->addChild(lambda_h);
    auto n_m_b_pole = std::make_shared<OperatorNode>("m_b_pole", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.m_b_pole = values[0]; return cache.m_b_pole; });
    n_m_b_pole->addChild(m_b_pole);
    auto n_m_c_pole = std::make_shared<OperatorNode>("m_c_pole", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.m_c_pole = values[0]; return cache.m_c_pole; });
    n_m_c_pole->addChild(m_c_pole);
    auto alpha_s_mu_b = std::make_shared<OperatorNode>("alpha_s(mu_b)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(w_config.hadronic_scale, MassType::POLE, MassType::POLE)); return cache.alpha_s_mu_b; });
    auto alpha_s_1_GeV = std::make_shared<OperatorNode>("alpha_s(1 GeV)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.alpha_s_1_GeV = ObsQCDProxy()(AlphasConfig(1.0, MassType::POLE, MassType::POLE)); return cache.alpha_s_1_GeV; });
    auto alpha_s_mu_f = std::make_shared<OperatorNode>("alpha_s(mu_f)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.alpha_s_mu_f = ObsQCDProxy()(AlphasConfig(values[0], MassType::POLE, MassType::POLE)); return cache.alpha_s_mu_f; });
    alpha_s_mu_f->addChild(mu_f);
    auto eta_f = std::make_shared<OperatorNode>("eta_f", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.eta_f = values[0] / values[1]; return cache.eta_f; });
    eta_f->addChildren({alpha_s_mu_f, alpha_s_1_GeV});
    auto alpha_s_mb_pole = std::make_shared<OperatorNode>("alpha_s(mb_pole)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.alpha_s_mb_pole = ObsQCDProxy()(AlphasConfig(cache.m_b_pole, MassType::POLE, MassType::POLE)); return cache.alpha_s_mb_pole; });
    alpha_s_mb_pole->addChild(n_m_b_pole);
    auto m_b_PS = std::make_shared<OperatorNode>("m_b_PS", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.m_b_PS = cache.m_b_pole - 4 * cache.alpha_s_mb_pole * cache.mu_f / (3 * PI); return cache.m_b_PS; });
    m_b_PS->addChildren({n_m_b_pole, alpha_s_mb_pole, mu_f});
    auto m_b_mu_b = std::make_shared<OperatorNode>("m_b(mu_b)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.m_b_mu_b = ObsQCDProxy()(MassConfig(5, w_config.hadronic_scale, MassType::POLE, MassType::POLE)); return cache.m_b_mu_b; });
    auto zc = std::make_shared<OperatorNode>("z_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.z_c = pow(cache.m_c_pole / cache.m_b_PS, 2); return cache.z_c; });
    zc->addChildren({n_m_c_pole, m_b_PS});
    auto Lb = std::make_shared<OperatorNode>("L_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.L_b = log(w_config.hadronic_scale / values[0]); return cache.L_b; });
    Lb->addChildren({m_b_PS});
    auto delta_M = std::make_shared<OperatorNode>("Delta_M", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.Delta_M = -6. * values[0] - 4. * (1. - values[1] / values[2]); return cache.Delta_M; });
    delta_M->addChildren({Lb, mu_f, m_b_PS});
    auto omega_0 = std::make_shared<OperatorNode>("omega_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.omega_0 = 2. * (values[0] - values[1]) / 3.; return cache.omega_0; });
    omega_0->addChildren({m_B, m_b_PS});
    auto lambda_u_hat = std::make_shared<OperatorNode>("lambda_u_hat", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.lambda_hat_u = values[0] * std::conj(values[1]) / values[2] * std::conj(values[3]); return cache.lambda_hat_u; });
    lambda_u_hat->addChildren({V_ub, V_us, V_tb, V_ts});
    auto f_K_perp = std::make_shared<OperatorNode>("f_K_perp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.f_K_perp = values[0] * std::pow(cache.eta_f, cache.C_F / cache.beta_0); return cache.f_K_perp; });
    f_K_perp->addChildren({f_K_perp_1GeV, eta_f, C_F, beta_0});
    auto a_1_perp = std::make_shared<OperatorNode>("a_1_perp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.a_1_perp = values[0] * std::pow(cache.eta_f, 4 * cache.C_F * (psi(2) + GAMMA - 0.5) / cache.beta_0); return cache.a_1_perp; });
    a_1_perp->addChildren({a_1_perp_1GeV, eta_f, C_F, beta_0});
    auto a_2_perp = std::make_shared<OperatorNode>("a_2_perp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.a_2_perp = values[0] * std::pow(cache.eta_f, 4 * cache.C_F * (psi(3) + GAMMA - 0.6667) / cache.beta_0); return cache.a_2_perp; });
    a_2_perp->addChildren({a_2_perp_1GeV, eta_f, C_F, beta_0});
    auto a_1_par = std::make_shared<OperatorNode>("a_1_par", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.a_1_par = values[0] * std::pow(cache.eta_f, 4 * cache.C_F * (psi(3) + GAMMA - 0.8333) / cache.beta_0); return cache.a_1_par; });
    a_1_par->addChildren({a_1_perp_1GeV, eta_f, C_F, beta_0});
    auto a_2_par = std::make_shared<OperatorNode>("a_2_par", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.a_2_par = values[0] * std::pow(cache.eta_f, 4 * cache.C_F * (psi(4) + GAMMA - 0.7917) / cache.beta_0); return cache.a_2_par; });
    a_2_par->addChildren({a_2_perp_1GeV, eta_f, C_F, beta_0});
    auto lambda_B_p = std::make_shared<OperatorNode>("lambda_B_p", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.lambda_B_p = values[0] / (1 - cache.alpha_s_mu_f * log(pow(cache.mu_f, 2)) * 1.8 / (3 * PI)); return cache.lambda_B_p; });
    lambda_B_p->addChildren({lambda_B_1GeV, alpha_s_mu_f});
    auto N_perp = std::make_shared<OperatorNode>("N_perp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.pref_T_perp = PI2 * values[0] * cache.f_K_perp / (cache.Nc * values[1]); return cache.pref_T_perp; });
    N_perp->addChildren({f_B, m_B, N_c, f_K_perp});
    auto N_par = std::make_shared<OperatorNode>("N_par", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.pref_T_par = PI2 * values[0] * values[1] * values[2] / (cache.Nc * values[3]); return cache.pref_T_par; });
    N_par->addChildren({f_B, f_K_par, m_Ks, m_B});
    auto T_par_m_0 = std::make_shared<OperatorNode>("T_par_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        double delta_qu = cfg.charge == Charge::B_PLUS ? 1. : 0.;
        cache.T_par_m_0 = 4 * values[0] / cache.m_b_PS * (3. * cache.lambda_hat_u * delta_qu * cache.C[WCoef::C2] - cache.C_bar[WCoef::C3] - 3. * cache.C_bar[WCoef::C4]);
        return cache.T_par_m_0; 
    });
    T_par_m_0->addChildren({m_B, m_b_PS, lambda_u_hat, wilson_cache, wilson_bar_cache});
    auto N_0 = std::make_shared<OperatorNode>("N_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.N_0 = values[0] * std::conj(values[1]) * values[2] / (values[3] * std::sqrt(3072. * std::pow(PI, 5) * std::pow(values[4], 3))); return cache.N_0; });
    N_0->addChildren({V_tb, V_ts, G_F, inv_alpha_em, m_B});

    auto T_perp_p_cache = std::make_shared<OperatorNode>("T_perp_p_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_perp_p, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8, std::placeholders::_9, std::placeholders::_10);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_perp_p_lookup, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]); 
        return T_perp_p_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_perp_p_cache->addChildren({m_B, m_Ks, f_B, f_K_par, zeta_3_A, zeta_3_V, omega_10_A, delta_tilde_p, delta_tilde_m, q2_min, n_q2_high, alpha_s_mu_b, N_perp, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_perp, a_2_perp, eq, lambda_B_p, wilson_cache, wilson_bar_cache});

    auto T_perp_m_cache = std::make_shared<OperatorNode>("T_perp_m_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_perp_m, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8, std::placeholders::_9, std::placeholders::_10);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_perp_m_lookup, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]); 
        return T_perp_m_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_perp_m_cache->addChildren({m_B, m_Ks, f_B, f_K_par, zeta_3_A, zeta_3_V, omega_10_A, delta_tilde_p, delta_tilde_m, q2_min, n_q2_high, alpha_s_mu_b, N_perp, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_perp, a_2_perp, eq, lambda_B_p, wilson_cache, wilson_bar_cache});

    auto T_par_p_cache = std::make_shared<OperatorNode>("T_par_p_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_par_p, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_par_p_lookup, values[0], values[1]); 
        return T_par_p_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_par_p_cache->addChildren({m_B, m_Ks, q2_min, n_q2_high, alpha_s_mu_b, N_par, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_par, a_2_par, eq, lambda_B_p, omega_0, T_par_m_0, wilson_cache, wilson_bar_cache});

    auto T_par_m_cache = std::make_shared<OperatorNode>("T_par_m_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_par_m, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_par_m_lookup, values[0], values[1]); 
        return T_par_m_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_par_m_cache->addChildren({m_B, m_Ks, q2_min, n_q2_high, alpha_s_mu_b, N_par, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_par, a_2_par, eq, lambda_B_p, omega_0, T_par_m_0, wilson_cache, wilson_bar_cache});

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

        double q2_min = 4 * std::pow(0.1057, 2);
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

    auto test_int = std::make_shared<OperatorNode>("test_int", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {
        double m_B = values[0];
        double m_K = values[1];
        double q2 = 0.99 * cache.q2_max;
        auto f = [m_B, m_K, q2, this] (double u) {
            double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
            double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
            complex_t i1 = phi * (T_par_p_p_f(u, q2, m_B, m_K) + T_par_p_nf(u, q2, m_B, m_K));
            complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K));
            return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2, m_B) * i2;
        };

        std::ofstream fs;
        fs.open("B_Ksll_int.csv");
        fs << "q2,f_re,f_im\n";

        auto write_line = [&] (double u) {
            fs << u 
            << "," << std::real(f(u)) << "," << std::imag(f(u))
            << "\n";
        };

        size_t n = 1000;
        double u0 {0}, u1 {1 - 1e-6};
        double du = (u1 - u0) / (n - 1);
        double u = u0;
        for (size_t i = 0; i < n; i++) {
            write_line(u);
            u += du;
        }

        return 0; 
    });
    test_int->addChildren({m_B, m_Ks, q2_min, q2_max, alpha_s_mu_b, N_par, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_par, a_2_par, eq, lambda_B_p, omega_0, T_par_m_0, wilson_cache, wilson_bar_cache});

    auto test_T = std::make_shared<OperatorNode>("test_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {
        std::ofstream fs;
        fs.open("B_Ksll_T.csv");
        fs << "q2,T_perp_p_re,T_perp_p_im,T_perp_m_re,T_perp_m_im,T_par_p_re,T_par_p_im,T_par_m_re,T_par_m_im\n";

        auto write_line = [&] (double q2) {
            fs << q2 
            << "," << std::real(T_perp_p_cached(q2)) << "," << std::imag(T_perp_p_cached(q2))
            << "," << std::real(T_perp_m_cached(q2)) << "," << std::imag(T_perp_m_cached(q2))
            << "," << std::real(T_par_p_cached(q2)) << "," << std::imag(T_par_p_cached(q2))
            << "," << std::real(T_par_m_cached(q2)) << "," << std::imag(T_par_m_cached(q2))
            << "\n";
        };

        size_t n = 200;
        double dq2 = (values[0] - cache.q2_min) / n;
        double q2 = cache.q2_min;
        for (size_t i = 0; i <= n; i++) {
            write_line(q2);
            q2 += dq2;
        }

        return 0; 
    });
    test_T->addChildren({q2_high, test_ff, T_perp_p_cache, T_perp_m_cache, T_par_p_cache, T_par_m_cache});

    auto test_J = std::make_shared<OperatorNode>("test_J", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {
        double m_B = values[0];
        double m_K = values[1];
        double m_l = values[2];
        double m_s = values[3];
        double f_B = values[4];
        double f_K_par = values[5];
        
        std::ofstream fs;
        fs.open("B_Ksll_J.csv");
        fs << "q2,J1s,J1c,J2s,J2c,J3,J4,J5,J6s,J6c,J7,J8,J9\n";

        auto write_line = [&] (double q2) {
            fs << q2 
            << "," << J1s(q2, m_B, m_K, m_l)
            << "," << J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)
            << "," << J2s(q2, m_B, m_K, m_l)
            << "," << J2c(q2, m_B, m_K, m_l)
            << "," << J3(q2, m_B, m_K, m_l)
            << "," << J4(q2, m_B, m_K, m_l)
            << "," << J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)
            << "," << J6s(q2, m_B, m_K, m_l)
            << "," << J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)
            << "," << J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s)
            << "," << J8(q2, m_B, m_K, m_l)
            << "," << J9(q2, m_B, m_K, m_l)
            << "\n";
        };

        size_t n = 200;
        double dq2 = (values[6] - cache.q2_min) / n;
        double q2 = cache.q2_min;
        for (size_t i = 0; i <= n; i++) {
            write_line(q2);
            q2 += dq2;
        }

        return 0; 
    });
    test_J->addChildren({m_B, m_Ks, m_l, m_s, f_B, f_K_par, q2_high, m_b_PS, m_b_mu_b, alpha_s_mu_b, C_F, N_c, lambda_B_p, test_ff, T_perp_p_cache, T_perp_m_cache, T_par_m_cache});

    roots.emplace(Observables::TEST_B__KS_L_L, test_J);
}