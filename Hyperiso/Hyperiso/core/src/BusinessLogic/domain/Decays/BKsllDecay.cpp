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
    return 1. / (1 - q2 / std::pow(m_R, 2));
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
    return (cache.tp * (mB2 - mK2 - q2) * A_1 - 16. * m_B * mK2 * (m_B + m_K) * A_12) / ((cache.tp - q2) * (cache.tm - q2));
}

double BKstarllDecay::T_3(double q2, double m_B, double m_K) {
    double T_2 = F_a(FF::T2, q2);
    double T_23 = F_a(FF::T23, q2);
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    return ((mB2 - mK2) * (mB2 + 3. * mK2 - q2) * T_2 - 8. * m_B * mK2 * (m_B - m_K) * T_23) / ((cache.tp - q2) * (cache.tm - q2));
}

double BKstarllDecay::xi_perp(double q2, double m_B, double m_K) {
    return m_B * F_a(FF::V, q2) / (m_B + m_K);
}

double BKstarllDecay::xi_par(double q2, double m_B, double m_K) {
    return (m_B + m_K) * F_a(FF::A1, q2) / (2. * E_K(q2, m_B, m_K)) - (m_B - m_K) * A_2(q2, m_B, m_K) / m_B;
}

double BKstarllDecay::f_perp(double q2, double m_B, double m_K) {
    return std::sqrt(2. * lambda(q2, m_B, m_K)) / (m_B + m_K) * F_a(FF::V, q2);
}

double BKstarllDecay::f_par(double q2, double m_B, double m_K) {
    return RT2 * (m_B + m_K) * F_a(FF::A1, q2);
}

double BKstarllDecay::f_0(double q2, double m_B, double m_K) {
    return ((m_B * m_B - q2 - m_K * m_K) * std::pow(m_B + m_K, 2) * F_a(FF::A1, q2) - lambda(q2, m_B, m_K) * A_2(q2, m_B, m_K)) / (2. * m_K * (m_B + m_K) * sqrt(q2));
}

double BKstarllDecay::F_perp(double s) {
    if (fpeq(s, 0.0)) 
        return 1.0 + cache.a_1_perp + 2.0 * cache.a_2_perp;
    
    double d = s - 1.;
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double s2 = s * s;
    double ls = std::log(s);
    double f0 = (s + 1.) / d2 - 2. * s * ls / d3;
    double f1 = -(s2 + 10. * s + 1.) / d3 + 6. * s * (s + 1) * ls / d4;
    double f2 = (s + 1.) + (s2 + 28. * s + 1.) / d4 - 12. * s * (s2 + 3. * s + 1.) * ls / d5;
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
    double f0 = (s - 3.) / d2 + 2. * ls / d3;
    double f1 = -(s2 - 8. * s - 17.) / d3 + 6. * (3. * s + 1.) * ls / d4;
    double f2 = -(s * s2 - 15. * s2 - 123. * s - 43.) / d4 - 12. * (6. * s2 + 8. * s + 1.) * ls / d5;
    return f0 + cache.a_1_perp * f1 + cache.a_2_perp * f2;
}

double BKstarllDecay::gv_dga_4(double u, double z3a, double z3v, double w10a, double dtp, double dtm) {
    double a1 = -60. * z3a * (w10a + 4.) + 1680. * z3v;
    double a2 = 30. * z3a * (15. * w10a + 32.) - 12600. * z3v + 36. * cache.a_1_par - 72. * cache.a_2_par - 12.;
    double a3 = -100. * z3a * (9. * w10a + 8.) + 25200. * z3v - 48. * cache.a_1_par + 240. * cache.a_2_par;
    double a4 = 525. * z3a * w10a - 14700. * z3v - 180. * cache.a_2_par;
    return -u * (a1 + u * (a2 + u * (a3 + u * a4))) / 4. + dtp * (9. * u - 1.5) + dtm * 6. * u + 3. * (dtp + dtm) * log(1 - u);
}

complex_t BKstarllDecay::F_V(double v, bool bar) {
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    return .75 * (
        h(v, cache.m_c_pole, w_config.hadronic_scale) * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6] + l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
      + h(v, cache.m_b_pole, w_config.hadronic_scale) * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6]) 
      + h(v, 0, w_config.hadronic_scale) * (cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4] + 3. * cache.C_bar[WCoef::C6] - l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
      - 8. / 27. * (cache.C_bar[WCoef::C3] - cache.C_bar[WCoef::C5] - 15. * cache.C_bar[WCoef::C6])
    );
}

double BKstarllDecay::L(double q2) {
    double mb2 = std::pow(cache.m_b_PS, 2);
    return (q2 - mb2) * std::log(1 - q2 / mb2) / q2;
}

complex_t BKstarllDecay::C_perp_0(double q2, double m_B, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return C7 + q2 * (Y(q2) + cache.lambda_hat_u * Y_u(q2)) / (2. * cache.m_b_PS * m_B);
}

complex_t BKstarllDecay::C_par_0(double q2, double m_B, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 - m_B * (Y(q2) + cache.lambda_hat_u * Y_u(q2)) / (2. * cache.m_b_PS);
}

complex_t BKstarllDecay::C_perp_f(double q2, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return C7 * (2. * std::log(cache.m_b_PS / w_config.hadronic_scale) - L(q2) + cache.Delta_M);
}

complex_t BKstarllDecay::C_par_f(double q2, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 * (2. * std::log(cache.m_b_PS / w_config.hadronic_scale) + 2. * L(q2) + cache.Delta_M);
}

complex_t BKstarllDecay::C_perp_nf(double q2, double m_B, bool bar) {
    double s_hat = q2 / (cache.m_b_PS * cache.m_b_PS);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t F_27 = f_27(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_27_u(s_hat, cache.L_b) * l_u;
    complex_t F_19 = f_19_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_19_u(s_hat, cache.L_b) * l_u;
    complex_t F_29 = f_29_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_29_u(s_hat, cache.L_b) * l_u;
    return -(
        cache.C_bar[WCoef::C2] * F_27 
      + cache.C[WCoef::C8] * f_87(s_hat, cache.L_b)
      + q2 / (2. * cache.m_b_PS * m_B) * (
            (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C1] / 3.) * F_29
          + 2. * cache.C_bar[WCoef::C1] * F_19
          + cache.C[WCoef::C8] * f_89(s_hat)
        )
    ) / cache.C_F;
}

complex_t BKstarllDecay::C_par_nf(double q2, double m_B, bool bar) {
    double s_hat = q2 / (m_B * m_B);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t F_27 = f_27(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_27_u(s_hat, cache.L_b) * l_u;
    complex_t F_19 = f_19_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_19_u(s_hat, cache.L_b) * l_u;
    complex_t F_29 = f_29_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_29_u(s_hat, cache.L_b) * l_u;
    return (
        cache.C_bar[WCoef::C2] * F_27
      + cache.C[WCoef::C8] * f_87(s_hat, cache.L_b)
      + m_B / (2 * cache.m_b_PS) * (
            (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C1] / 3.) * F_29
           + 2. * cache.C_bar[WCoef::C1] * F_19
           + cache.C[WCoef::C8] * f_89(s_hat))
    ) / cache.C_F;
}

complex_t BKstarllDecay::T_par_m_0(double m_B, bool bar) {
    int delta_qu = cfg.charge == Charge::B_PLUS;
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    return 4. * m_B / cache.m_b_PS * (delta_qu * 3. * l_u * cache.C[WCoef::C2] - cache.C_bar[WCoef::C3] - 3. * cache.C_bar[WCoef::C4]);
}

complex_t BKstarllDecay::T_par_p_p_f(double u, double q2, double m_B, double m_K, bool bar) {
    return 2. * T_perp_p_p_f(u, q2, m_B, m_K, bar);
}

complex_t BKstarllDecay::T_par_p_m_f(double u, double q2, double m_B, double m_K, bool bar) {
    return 2. * T_perp_p_m_f(u, q2, m_B, m_K, bar);
}

complex_t BKstarllDecay::T_perp_p_p_f(double u, double q2, double m_B, double m_K, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return 2. * m_B / (1. - u) / E_K(q2, m_B, m_K) * C7;
}

complex_t BKstarllDecay::T_perp_p_m_f(double u, double q2, double m_B, double m_K, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] - cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return 2. * m_B / (1. - u) / E_K(q2, m_B, m_K) * C7;
} 

complex_t BKstarllDecay::T_perp_p_nf(double u, double q2, double m_B, double m_K, bool bar) {
    double E = E_K(q2, m_B, m_K);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t t_perp_mc = t_perp(u, cache.m_c_pole, q2, E, m_B);
    complex_t t_perp_mb = t_perp(u, cache.m_b_PS, q2, E, m_B);
    complex_t t_perp_0 = t_perp(u, 0, q2, E, m_B);
    return -4 * cache.e_d * cache.C[WCoef::C8] / (u + (1 - u) * q2 / (m_B * m_B))
            + m_B / (2 * cache.m_b_PS) * (
                cache.e_u * (
                    t_perp_mc * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] + l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.))
                    - t_perp_0 * l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)
                )
                + cache.e_d * (
                    t_perp_mb * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] - 4 * cache.m_b_PS / m_B * cache.C_bar[WCoef::C5])
                    + t_perp_0 * cache.C_bar[WCoef::C3]
                )
            );
}

complex_t BKstarllDecay::T_par_p_nf(double u, double q2, double m_B, double m_K, bool bar) {
    double E = E_K(q2, m_B, m_K);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t t_par_mc = t_par(u, cache.m_c_pole, q2, E, m_B);
    complex_t t_par_mb = t_par(u, cache.m_b_PS, q2, E, m_B);
    complex_t t_par_0 = t_par(u, 0, q2, E, m_B);
    return m_B / cache.m_b_PS * (
        cache.e_u * (
            t_par_mc * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] + l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.))
            - t_par_0 * l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)
        )
        + cache.e_d * (
            t_par_mb * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6])
            + t_par_0 * cache.C_bar[WCoef::C3]
        )
    );
}

complex_t BKstarllDecay::T_par_m_nf(double u, double q2, double m_B, double m_K, bool bar) {
    double v = m_B * m_B * (1 - u) + q2 * u;
    return 8. * m_B * m_B * cache.C[WCoef::C8] / v + 8. * m_B / cache.m_b_PS * F_V(v, bar);
}

complex_t BKstarllDecay::inv_lambda_B_m(double q2, double m_B) {
    double x = q2 / (m_B * cache.omega_0);
    return std::exp(-x) / cache.omega_0 * (I * PI - Ei(x));
}

complex_t BKstarllDecay::I_perp_p(double q2, double m_B, double m_K, bool bar) {
    auto f = [q2, m_B, m_K, bar, this] (double u) {
        return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * (T_perp_p_p_f(u, q2, m_B, m_K, bar) + T_perp_p_nf(u, q2, m_B, m_K, bar));
    };
    
    return cache.alpha_s_mu_f / (4. * PI) * cache.C_F / cache.lambda_B_p * c_integrate(f, 0, 1, 1e-2);
}

complex_t BKstarllDecay::I_perp_m(double q2, double m_B, double m_K, bool bar) {
    auto f = [q2, m_B, m_K, bar, this] (double u) {
        return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * (T_perp_p_m_f(u, q2, m_B, m_K, bar) + T_perp_p_nf(u, q2, m_B, m_K, bar));
    };
    
    return cache.alpha_s_mu_f / (4. * PI) * cache.C_F / cache.lambda_B_p * c_integrate(f, 0, 1, 1e-2);
}

complex_t BKstarllDecay::I_par_p(double q2, double m_B, double m_K, bool bar) {
    auto f = [q2, m_B, m_K, bar, this] (double u) {
        double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        complex_t i1 = phi * (T_par_p_p_f(u, q2, m_B, m_K, bar) + T_par_p_nf(u, q2, m_B, m_K, bar));
        complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K, bar));
        return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2, m_B) * i2;
    };
    
    return c_integrate(f, 0, 1 - 1e-6, 1e-2);
}

complex_t BKstarllDecay::I_par_m(double q2, double m_B, double m_K, bool bar) {
    auto f = [q2, m_B, m_K, bar, this] (double u) {
        double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        complex_t i1 = phi * (T_par_p_m_f(u, q2, m_B, m_K, bar) + T_par_p_nf(u, q2, m_B, m_K, bar));
        complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K, bar));
        return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2, m_B) * i2;
    };
    
    return c_integrate(f, 0, 1 - 1e-6, 1e-2);
}

complex_t BKstarllDecay::I_HSA_1(double q2, double m_B, bool bar) {
    auto f = [q2, m_B, bar, this] (double u) {
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        double v = m_B * m_B * (1 - u) + u * q2;
        return phi * m_B * m_B / v * F_V(v, bar);
    };
    
    return c_integrate(f, 0, 1, 1e-2);
}

complex_t BKstarllDecay::I_HSA_2(double q2, double m_B, double z3a, double z3v, double w10a, double dtp, double dtm, bool bar) {
    auto f = [q2, m_B, z3a, z3v, w10a, dtp, dtm, bar, this] (double u) {
        double int_phi_par = gv_dga_4(u, z3a, z3v, w10a, dtp, dtm);
        double v = m_B * m_B * (1 - u) + u * q2;
        return int_phi_par * F_V(v, bar);
    };
    
    return c_integrate(f, 0, 1, 1e-2);
}

complex_t BKstarllDecay::delta_T_perp_WA(double q2, double m_B, double m_K, double f_B, double f_K_par) {
    double pref = cache.e_q * 2. * PI2 * f_B / (cache.m_b_PS * m_B);
    complex_t W_perp = cache.C[WCoef::C3] + 4. / .3 * (cache.C[WCoef::C4] + 3. * cache.C[WCoef::C5] + 4. * cache.C[WCoef::C6]);
    complex_t W_par = cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4];
    if (cfg.charge == Charge::B_PLUS) 
        W_par += -3. * cache.C[WCoef::C2];
    double s_hat = q2 / (m_B * m_B);
    return pref * (
        -2. * cache.f_K_perp * W_perp * F_perp(s_hat)
      + f_K_par * m_K * W_par / (3. * (1 - s_hat) * cache.lambda_B_p)
    );
}

complex_t BKstarllDecay::delta_T_perp_HSA(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm, bool bar) {
    double pref = cache.e_q * cache.alpha_s_mu_b * cache.C_F * PI * f_B / (cache.Nc * cache.m_b_PS * m_B);
    double s_hat = q2 / (m_B * m_B);
    return pref * (
        3. * cache.C[WCoef::C8] * cache.m_b_PS / m_B * cache.f_K_perp * X_perp(s_hat)
      + 2. * cache.f_K_perp * I_HSA_1(q2, m_B, bar)
      - m_K * f_K_par / ((1 - s_hat) * cache.lambda_B_p) * I_HSA_2(q2, m_B, z3a, z3v, w10a, dtp, dtm, bar)
    );
}

complex_t BKstarllDecay::T_perp_p(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm, bool bar) {
    complex_t C_perp_p = C_perp_0(q2, m_B, 1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_perp_f(q2, 1, bar) + C_perp_nf(q2, m_B, bar));
    return xi_perp(q2, m_B, m_K) * C_perp_p + cache.pref_T_perp * I_perp_p(q2, m_B, m_K, bar) + delta_T_perp_WA(q2, m_B, m_K, f_B, f_K_par) + delta_T_perp_HSA(q2, m_B, m_K, f_B, f_K_par, z3a, z3v, w10a, dtp, dtm, bar);
}

complex_t BKstarllDecay::T_perp_m(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm, bool bar) {
    complex_t C_perp_m = C_perp_0(q2, m_B, -1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_perp_f(q2, -1, bar) + C_perp_nf(q2, m_B, bar));
    return xi_perp(q2, m_B, m_K) * C_perp_m + cache.pref_T_perp * I_perp_m(q2, m_B, m_K, bar) + delta_T_perp_WA(q2, m_B, m_K, f_B, f_K_par) + delta_T_perp_HSA(q2, m_B, m_K, f_B, f_K_par, z3a, z3v, w10a, dtp, dtm, bar);
}

complex_t BKstarllDecay::T_par_p(double q2, double m_B, double m_K, bool bar) {
    complex_t C_par_p = C_par_0(q2, m_B, 1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_par_f(q2, 1, bar) + C_par_nf(q2, m_B, bar));
    return xi_par(q2, m_B, m_K) * C_par_p + cache.pref_T_par / E_K(q2, m_B, m_K) * I_par_p(q2, m_B, m_K, bar);
}

complex_t BKstarllDecay::T_par_m(double q2, double m_B, double m_K, bool bar) {
    complex_t C_par_m = C_par_0(q2, m_B, -1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_par_f(q2, -1, bar) + C_par_nf(q2, m_B, bar));
    return xi_par(q2, m_B, m_K) * C_par_m + cache.pref_T_par / E_K(q2, m_B, m_K) * I_par_m(q2, m_B, m_K, bar);
}

complex_t BKstarllDecay::Delta_par(double q2, double m_B, double m_K, double f_B, double f_K_par) {
    return 1. + cache.alpha_s_mu_b * cache.C_F / (2. * PI) * (
        L(q2) - 1 - 3. * PI2 * q2 * f_B * f_K_par * m_K / (cache.Nc * m_B * cache.lambda_B_p * xi_par(q2, m_B, m_K) * std::pow(E_K(q2, m_B, m_K), 3)) * F_perp(0.0)
    );
}

complex_t BKstarllDecay::T_perp_p_cached(double q2, bool bar) {
    return lerp(q2, bar ? T_perp_p_bar_lookup : T_perp_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_perp_m_cached(double q2, bool bar) {
    return lerp(q2, bar ? T_perp_m_bar_lookup : T_perp_m_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_par_p_cached(double q2, bool bar) {
    return lerp(q2, bar ? T_par_p_bar_lookup : T_par_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_par_m_cached(double q2, bool bar) {
    return lerp(q2, bar ? T_par_m_bar_lookup : T_par_m_lookup, cache.q2_min, cache.q2_high);
}

double BKstarllDecay::beta_l(double q2, double m_l) {
    return std::sqrt(1 - 4. * m_l * m_l / q2);
}

double BKstarllDecay::lambda(double q2, double m_B, double m_K) {
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    return mB2 * mB2 + mK2 * mK2 + q2 * q2 - 2. * (mB2 * mK2 + (mB2 + mK2) * q2);
}

complex_t BKstarllDecay::N(double q2, double m_B, double m_K, double m_l, bool bar) {
    complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0; 
    return N0 * std::sqrt(q2 * beta_l(q2, m_l) * std::sqrt(lambda(q2, m_B, m_K)));
}

complex_t BKstarllDecay::A_perp_low(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
    if (bar) w = std::conj(w);
    return N(q2, m_B, m_K, m_l, bar) * std::sqrt(2 * lambda(q2, m_B, m_K)) * (
        w * F_a(FF::V, q2) / (m_B + m_K)
      + 2. * cache.m_b_PS * T_perp_p_cached(q2, bar) / q2
    );
}

complex_t BKstarllDecay::A_par_low(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
    if (bar) w = std::conj(w);
    return -N(q2, m_B, m_K, m_l, bar) * std::sqrt(2.) * (m_B * m_B - m_K * m_K) * (
        w * F_a(FF::A1, q2) / (m_B - m_K)
      + 4. * cache.m_b_PS * E_K(q2, m_B, m_K) * T_perp_m_cached(q2, bar) / (m_B * q2)
    );
}

complex_t BKstarllDecay::A_0_low(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    double mB2 = m_B * m_B;
    double mK2 = m_K * m_K;
    complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
    if (bar) w = std::conj(w);
    return -N(q2, m_B, m_K, m_l, bar) / (2. * m_K * std::sqrt(q2)) * (
        w * (
            (mB2 - mK2 - q2) * (m_B + m_K) * F_a(FF::A1, q2)
          - lambda(q2, m_B, m_K) * A_2(q2, m_B, m_K) / (m_B + m_K)
        )
      + 2. * cache.m_b_PS * (
          2. * E_K(q2, m_B, m_K) * (mB2 + 3. * mK2 - q2) * T_perp_m_cached(q2, bar) 
        - lambda(q2, m_B, m_K) * (T_perp_m_cached(q2, bar) + T_par_m_cached(q2, bar)) / (mB2 - mK2)
      )
    );
}

complex_t BKstarllDecay::A_t_low(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }
    return N(q2, m_B, m_K, m_l, bar) * std::sqrt(lambda(q2, m_B, m_K) / q2) * (
        (C10 + q2 / (m_l * (cache.m_b_mu_b + m_s)) * CQ2)
    ) * E_K(q2, m_B, m_K) * xi_par(q2, m_B, m_K) / (m_K * Delta_par(q2, m_B, m_K, f_B, f_K_par));
}

complex_t BKstarllDecay::A_S_low(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);
    return -2. * N(q2, m_B, m_K, m_l, bar) / (cache.m_b_mu_b + m_s) * std::sqrt(lambda(q2, m_B, m_K)) 
            * CQ1 * E_K(q2, m_B, m_K) * xi_par(q2, m_B, m_K) / (m_K * Delta_par(q2, m_B, m_K, f_B, f_K_par));
}

complex_t BKstarllDecay::C7_eff(double q2, double m_B, bool bar) {
    complex_t A = A_Seidel(q2, cache.m_b_PS, w_config.hadronic_scale);
    return (bar ? std::conj(cache.C[WCoef::C7]) : cache.C[WCoef::C7]) + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*f_87(q2 / (m_B * m_B), cache.L_b));
}

complex_t BKstarllDecay::C9_eff(double q2, double m_B, bool bar) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+l_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    complex_t A = A_Seidel(q2, cache.m_b_PS, w_config.hadronic_scale);
    complex_t B = B_Seidel(q2, cache.m_b_PS, w_config.hadronic_scale);
    complex_t C = C_Seidel(q2, w_config.hadronic_scale);

    return (bar ? std::conj(cache.C[WCoef::C9]) : cache.C[WCoef::C9])
         + h(q2, 0., w_config.hadronic_scale) * C_h0
         + h(q2, cache.m_b_PS, w_config.hadronic_scale) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * f_89(q2 / (m_B * m_B)))
         + std::pow(cache.m_c_mu_b, 2) / q2 * C_mc;
}

complex_t BKstarllDecay::A_perp_high(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, m_B, bar) + (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, m_B, bar) + (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] + cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return N(q2, m_B, m_K, m_l, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * m_B / q2 * C7) * f_perp(q2, m_B, m_K);
}

complex_t BKstarllDecay::A_par_high(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, m_B, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, m_B, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, m_B, m_K, m_l, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * m_B / q2 * C7) * f_par(q2, m_B, m_K);
}

complex_t BKstarllDecay::A_0_high(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, m_B, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, m_B, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, m_B, m_K, m_l, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * m_B / q2 * C7) * f_0(q2, m_B, m_K);
}

complex_t BKstarllDecay::A_t_high(double q2, double m_B, double m_K, double m_l, double m_s, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }
    return N(q2, m_B, m_K, m_l, bar) / sqrt(q2 * lambda(q2, m_B, m_K)) * (2. * C10 + q2 / m_l * CQ2 / (cache.m_b_mu_b + m_s)) * F_a(FF::A0, q2);
}

complex_t BKstarllDecay::A_S_high(double q2, double m_B, double m_K, double m_l, double m_s, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);
    return -2. * N(q2, m_B, m_K, m_l, bar) * sqrt(lambda(q2, m_B, m_K)) * CQ1 / (cache.m_b_mu_b + m_s) * F_a(FF::A0, q2);
}

complex_t BKstarllDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BKstarllDecay::A_perp(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    return interpolate(q2, A_perp_low(q2, m_B, m_K, m_l, sign, bar), A_perp_high(q2, m_B, m_K, m_l, sign, bar));
}

complex_t BKstarllDecay::A_par(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    return interpolate(q2, A_par_low(q2, m_B, m_K, m_l, sign, bar), A_par_high(q2, m_B, m_K, m_l, sign, bar));
}

complex_t BKstarllDecay::A_0(double q2, double m_B, double m_K, double m_l, double sign, bool bar) {
    return interpolate(q2, A_0_low(q2, m_B, m_K, m_l, sign, bar), A_0_high(q2, m_B, m_K, m_l, sign, bar));
}

complex_t BKstarllDecay::A_t(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    return interpolate(q2, A_t_low(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar), A_0_high(q2, m_B, m_K, m_l, m_s, bar));
}

complex_t BKstarllDecay::A_S(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    return interpolate(q2, A_S_low(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar), A_S_high(q2, m_B, m_K, m_l, m_s, bar));
}

double BKstarllDecay::J1s(double q2, double m_B, double m_K, double m_l, bool bar) {
    return (2. + std::pow(beta_l(q2, m_l), 2)) / 4. * (
        std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, bar)), 2)
    ) + std::pow(2. * m_l, 2) / q2 * std::real(
        A_perp(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar))
      + A_perp(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar))
    );
}

double BKstarllDecay::J1c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    return std::pow(std::abs(A_0(q2, m_B, m_K, m_l, -1, bar)), 2) + std::pow(std::abs(A_0(q2, m_B, m_K, m_l, 1, bar)), 2) 
         + std::pow(2 * m_l, 2) / q2 * (
              std::pow(std::abs(A_t(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar)), 2)
            + 2. * std::real(A_0(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_0(q2, m_B, m_K, m_l, 1, bar)))
           )
         + std::pow(beta_l(q2, m_l) * std::abs(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar)), 2);
}

double BKstarllDecay::J2s(double q2, double m_B, double m_K, double m_l, bool bar) {
    return std::pow(beta_l(q2, m_l), 2) / 4. * (
        std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, bar)), 2)
    );
}

double BKstarllDecay::J2c(double q2, double m_B, double m_K, double m_l, bool bar) {
    return -std::pow(beta_l(q2, m_l), 2) * (
        std::pow(std::abs(A_0(q2, m_B, m_K, m_l, -1, bar)), 2) 
      + std::pow(std::abs(A_0(q2, m_B, m_K, m_l, 1, bar)), 2)
    );
}

double BKstarllDecay::J3(double q2, double m_B, double m_K, double m_l, bool bar) {
    return std::pow(beta_l(q2, m_l), 2) / 2. * (
        std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, bar)), 2)
      - std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1, bar)), 2)
      - std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1, bar)), 2)
    );
}

double BKstarllDecay::J4(double q2, double m_B, double m_K, double m_l, bool bar) {
    return std::pow(beta_l(q2, m_l), 2) / std::sqrt(2.) * (
        std::real(A_0(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, bar))) 
      + std::real(A_0(q2, m_B, m_K, m_l, 1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar)))
    );
}

double BKstarllDecay::J5(double q2, double m_B,double m_K,double m_l,double f_B,double f_K_par,double m_s, bool bar) {
    return beta_l(q2, m_l) * std::sqrt(2.) * (
        std::real(A_0(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, bar))) 
      - std::real(A_0(q2, m_B, m_K, m_l, 1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar)))
      - m_l / std::sqrt(q2) * std::real((A_perp(q2, m_B, m_K, m_l, -1, bar) + A_perp(q2, m_B, m_K, m_l, 1, bar)) * std::conj(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar)))
    );
}

double BKstarllDecay::J6s(double q2, double m_B, double m_K, double m_l, bool bar) {
    return 2. * beta_l(q2, m_l) * (
        std::real(A_perp(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, bar))) 
      - std::real(A_perp(q2, m_B, m_K, m_l, 1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar))) 
    );
}

double BKstarllDecay::J6c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    return 4. * beta_l(q2, m_l) * m_l / std::sqrt(q2) * (std::real((A_0(q2, m_B, m_K, m_l, -1, bar) + A_0(q2, m_B, m_K, m_l, 1, bar)) * std::conj(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar))));
}

double BKstarllDecay::J7(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    return beta_l(q2, m_l) * std::sqrt(2.) * (
        std::imag(A_0(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, bar))) 
      - std::imag(A_0(q2, m_B, m_K, m_l, 1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar)))
      + m_l / std::sqrt(q2) * std::imag((A_perp(q2, m_B, m_K, m_l, -1, bar) + A_perp(q2, m_B, m_K, m_l, 1, bar)) * std::conj(A_S(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar)))
    );
}

double BKstarllDecay::J8(double q2, double m_B, double m_K, double m_l, bool bar) {
    return std::pow(beta_l(q2, m_l), 2) / std::sqrt(2.) * (
        std::imag(A_0(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, bar))) 
      + std::imag(A_0(q2, m_B, m_K, m_l, 1, bar) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, bar)))
    );
}

double BKstarllDecay::J9(double q2, double m_B, double m_K, double m_l, bool bar) {
    return std::pow(beta_l(q2, m_l), 2) * (
        std::imag(A_perp(q2, m_B, m_K, m_l, -1, bar) * std::conj(A_par(q2, m_B, m_K, m_l, -1, bar))) 
      + std::imag(A_perp(q2, m_B, m_K, m_l, 1, bar) * std::conj(A_par(q2, m_B, m_K, m_l, 1, bar))) 
    );
}

void BKstarllDecay::compute_binned_J_i(double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    auto fill_binned = [&] (std::array<std::vector<double>, 14>& dest, bool bar) {
        for (auto [q2_l, q2_u] : this->bins) {
            dest[0].emplace_back(integrate([&] (double q2) { return 2 * J1s(q2, m_B, m_K, m_l, bar) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar); }, q2_l, q2_u, 1e-2));
            dest[1].emplace_back(integrate([&] (double q2) { return J2s(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[2].emplace_back(integrate([&] (double q2) { return J2c(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[3].emplace_back(integrate([&] (double q2) { return J3(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[4].emplace_back(integrate([&] (double q2) { return J4(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[5].emplace_back(integrate([&] (double q2) { return J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar); }, q2_l, q2_u, 1e-2));
            dest[6].emplace_back(integrate([&] (double q2) { return beta_l(q2, m_l) * J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar); }, q2_l, q2_u, 1e-2));
            dest[7].emplace_back(integrate([&] (double q2) { return J6s(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[8].emplace_back(integrate([&] (double q2) { return beta_l(q2, m_l) * J6s(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[9].emplace_back(integrate([&] (double q2) { return J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar); }, q2_l, q2_u, 1e-2));
            dest[10].emplace_back(integrate([&] (double q2) { return J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar); }, q2_l, q2_u, 1e-2));
            dest[11].emplace_back(integrate([&] (double q2) { return beta_l(q2, m_l) * J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar); }, q2_l, q2_u, 1e-2));
            dest[12].emplace_back(integrate([&] (double q2) { return J8(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
            dest[13].emplace_back(integrate([&] (double q2) { return J9(q2, m_B, m_K, m_l, bar); }, q2_l, q2_u, 1e-2));
        }
    };

    fill_binned(this->J_i_binned, false);
    fill_binned(this->J_i_bar_binned, true);
}

double BKstarllDecay::dG_dq2_avg(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J1 = 2 * J1s(q2, m_B, m_K, m_l, false) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false);
    double J1bar = 2 * J1s(q2, m_B, m_K, m_l, true) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double J2 = 2 * J2s(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, false);
    double J2bar = 2 * J2s(q2, m_B, m_K, m_l, true) + J2c(q2, m_B, m_K, m_l, true);
    double dG = 0.75 * (J1 - J2 / 3.);
    double dGbar = 0.75 * (J1bar - J2bar / 3.);
    return dG + dGbar;
}

double BKstarllDecay::dG_dq2(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s, bool bar) {
    double J1 = 2 * J1s(q2, m_B, m_K, m_l, bar) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, bar);
    double J2 = 2 * J2s(q2, m_B, m_K, m_l, bar) + J2c(q2, m_B, m_K, m_l, bar);
    return 0.75 * (J1 - J2 / 3.);
}

double BKstarllDecay::A_FB(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J6 = 2 * J6s(q2, m_B, m_K, m_l, false) + J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false);
    double J6bar = 2 * J6s(q2, m_B, m_K, m_l, true) + J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return -0.375 * (J6 + J6bar) / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_CP(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J1 = 2 * J1s(q2, m_B, m_K, m_l, false) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false);
    double J1bar = 2 * J1s(q2, m_B, m_K, m_l, true) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double J2 = 2 * J2s(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, false);
    double J2bar = 2 * J2s(q2, m_B, m_K, m_l, true) + J2c(q2, m_B, m_K, m_l, true);
    double dG = 0.75 * (J1 - J2 / 3.);
    double dGbar = 0.75 * (J1bar - J2bar / 3.);
    return (dG - dGbar) / (dG + dGbar);
}

double BKstarllDecay::F_L(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J1 = J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double J2 = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    return 0.75 * (3 * J1 - J2) / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::F_T(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2 = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    return 4. * J2 / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_T_1(double q2, double m_B, double m_K, double m_l) {
    double AperpApar = std::real(A_par(q2, m_B, m_K, m_l, 1, false) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, false)) + A_par(q2, m_B, m_K, m_l, -1, false) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, false)));
    double AperpApar_bar = std::real(A_par(q2, m_B, m_K, m_l, 1, true) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, true)) + A_par(q2, m_B, m_K, m_l, -1, true) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, true)));
    double Aperp2 = std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, false)), 2) + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, false)), 2) + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, true)), 2) + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, true)), 2);
    double Apar2 = std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1, false)), 2) + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1, false)), 2) + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1, true)), 2) + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1, true)), 2);
    return -2. * (AperpApar + AperpApar_bar) / (Aperp2 + Apar2);
}

double BKstarllDecay::A_T_2(double q2, double m_B, double m_K, double m_l) {
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    return 0.5 * J3cp / J2scp;
}

double BKstarllDecay::A_T_3(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    double J4cp = J4(q2, m_B, m_K, m_l, false) + J4(q2, m_B, m_K, m_l, true);
    double J7cp = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double bl = beta_l(q2, m_l);
    return std::sqrt((4. * std::pow(J4cp, 2) + std::pow(bl * J7cp, 2)) / (-2. * J2ccp * (2. * J2scp + J3cp)));
}

double BKstarllDecay::A_T_4(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J4cp = J4(q2, m_B, m_K, m_l, false) + J4(q2, m_B, m_K, m_l, true);
    double J5cp = J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double J7cp = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double J8cp = J8(q2, m_B, m_K, m_l, false) + J8(q2, m_B, m_K, m_l, true);
    double bl = beta_l(q2, m_l);
    return std::sqrt((4. * std::pow(J8cp, 2) + std::pow(bl * J5cp, 2)) / (4. * std::pow(J4cp, 2) + std::pow(bl * J7cp, 2)));
}

double BKstarllDecay::A_T_5(double q2, double m_B, double m_K, double m_l) {
    complex_t AperpApar = A_perp(q2, m_B, m_K, m_l, 1, false) * std::conj(A_par(q2, m_B, m_K, m_l, -1, false)) + A_perp(q2, m_B, m_K, m_l, -1, false) * std::conj(A_par(q2, m_B, m_K, m_l, 1, false));
    complex_t AperpApar_bar = A_perp(q2, m_B, m_K, m_l, 1, true) * std::conj(A_par(q2, m_B, m_K, m_l, -1, true)) + A_perp(q2, m_B, m_K, m_l, -1, true) * std::conj(A_par(q2, m_B, m_K, m_l, 1, true));
    double Aperp2 = std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, false)), 2) + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, false)), 2) + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, 1, true)), 2) + std::pow(std::abs(A_perp(q2, m_B, m_K, m_l, -1, true)), 2);
    double Apar2 = std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1, false)), 2) + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1, false)), 2) + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, 1, true)), 2) + std::pow(std::abs(A_par(q2, m_B, m_K, m_l, -1, true)), 2);
    return std::abs(AperpApar + AperpApar_bar) / (Aperp2 + Apar2);
}

double BKstarllDecay::A_T_Re(double q2, double m_B, double m_K, double m_l) {
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J6scp = J6s(q2, m_B, m_K, m_l, false) + J6s(q2, m_B, m_K, m_l, true);
    double bl = beta_l(q2, m_l);
    return 0.25 * bl * J6scp / J2scp;
}

double BKstarllDecay::AA_T_Re(double q2, double m_B, double m_K, double m_l) {
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J6scpa = J6s(q2, m_B, m_K, m_l, false) - J6s(q2, m_B, m_K, m_l, true);
    double bl = beta_l(q2, m_l);
    return 0.25 * bl * J6scpa / J2scp;
}

double BKstarllDecay::A_Im(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J9cp = J9(q2, m_B, m_K, m_l, false) + J9(q2, m_B, m_K, m_l, true);
    return J9cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::alpha_K(double q2, double m_B, double m_K, double m_l) {
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J2cp = 2 * J2scp + J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    return -0.5 * J2cp / J2scp;
}

double BKstarllDecay::H_T_1(double q2, double m_B, double m_K, double m_l) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    double J4cp = J4(q2, m_B, m_K, m_l, false) + J4(q2, m_B, m_K, m_l, true);
    return RT2 * J4cp / std::sqrt(-J2ccp * (2 * J2scp - J3cp));
}

double BKstarllDecay::H_T_2(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    double J5cp = J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double bl = beta_l(q2, m_l);
    return bl * J5cp / std::sqrt(-2. * J2ccp * (2 * J2scp + J3cp));
}

double BKstarllDecay::H_T_3(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J6cp = 2 * (J6s(q2, m_B, m_K, m_l, false) + J6s(q2, m_B, m_K, m_l, true)) + J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    return 0.5 * J6cp / std::sqrt(4. * std::pow(J2scp, 2) - std::pow(J3cp, 2));
}

double BKstarllDecay::P_2(double q2, double m_B, double m_K, double m_l) {
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J6scp = J6s(q2, m_B, m_K, m_l, false) + J6s(q2, m_B, m_K, m_l, true);
    return 0.125 * J6scp / J2scp;
}

double BKstarllDecay::P_3(double q2, double m_B, double m_K, double m_l) {
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J9cp = J9(q2, m_B, m_K, m_l, false) + J9(q2, m_B, m_K, m_l, true);
    return -0.25 * J9cp / J2scp;
}

double BKstarllDecay::P_6(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    double J7cp = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    double bl = beta_l(q2, m_l);
    return -bl * J7cp / std::sqrt(-2. * J2ccp * (2 * J2scp - J3cp));
}

double BKstarllDecay::P_8(double q2, double m_B, double m_K, double m_l) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    double J8cp = J8(q2, m_B, m_K, m_l, false) + J8(q2, m_B, m_K, m_l, true);
    return -RT2 * J8cp / std::sqrt(-J2ccp * (2 * J2scp + J3cp));
}

double BKstarllDecay::Pp_4(double q2, double m_B, double m_K, double m_l) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J4cp = J4(q2, m_B, m_K, m_l, false) + J4(q2, m_B, m_K, m_l, true);
    return J4cp / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::Pp_5(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J5cp = J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return 0.5 * J5cp / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::Pp_6(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J7cp = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return -0.5 * J7cp / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::Pp_8(double q2, double m_B, double m_K, double m_l) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J8cp = J8(q2, m_B, m_K, m_l, false) + J8(q2, m_B, m_K, m_l, true);
    return -J8cp / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::S_3(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J3cp = J3(q2, m_B, m_K, m_l, false) + J3(q2, m_B, m_K, m_l, true);
    return J3cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::S_4(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J4cp = J4(q2, m_B, m_K, m_l, false) + J4(q2, m_B, m_K, m_l, true);
    return J4cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::S_5(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J5cp = J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return J5cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::S_6c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J6ccp = J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return J6ccp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::S_7(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J7cp = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) + J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return J7cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::S_8(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J8cp = J8(q2, m_B, m_K, m_l, false) + J8(q2, m_B, m_K, m_l, true);
    return J8cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::S_9(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J9cp = J9(q2, m_B, m_K, m_l, false) + J9(q2, m_B, m_K, m_l, true);
    return J9cp / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_3(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J3cpa = J3(q2, m_B, m_K, m_l, false) - J3(q2, m_B, m_K, m_l, true);
    return J3cpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_4(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J4cpa = J4(q2, m_B, m_K, m_l, false) - J4(q2, m_B, m_K, m_l, true);
    return J4cpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_5(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J5cpa = J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) - J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return J5cpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_6s(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J6scpa = J6s(q2, m_B, m_K, m_l, false) - J6s(q2, m_B, m_K, m_l, true);
    return J6scpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_7(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J7cpa = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) - J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return J7cpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_8(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J8cpa = J8(q2, m_B, m_K, m_l, false) - J8(q2, m_B, m_K, m_l, true);
    return J8cpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::A_9(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J9cpa = J9(q2, m_B, m_K, m_l, false) - J9(q2, m_B, m_K, m_l, true);
    return J9cpa / dG_dq2_avg(q2, m_B, m_K, m_l, f_B, f_K_par, m_s);
}

double BKstarllDecay::AP_1(double q2, double m_B, double m_K, double m_l) {
    double J3cpa = J3(q2, m_B, m_K, m_l, false) - J3(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    return 0.5 * J3cpa / J2scp;
}

double BKstarllDecay::AP_2(double q2, double m_B, double m_K, double m_l) {
    double J6scpa = J6s(q2, m_B, m_K, m_l, false) - J6s(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    return 0.125 * J6scpa / J2scp;
}

double BKstarllDecay::AP_3(double q2, double m_B, double m_K, double m_l) {
    double J9cpa = J9(q2, m_B, m_K, m_l, false) - J9(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    return -0.25 * J9cpa / J2scp;
}

double BKstarllDecay::APp_4(double q2, double m_B, double m_K, double m_l) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J4cpa = J4(q2, m_B, m_K, m_l, false) - J4(q2, m_B, m_K, m_l, true);
    return J4cpa / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::APp_5(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J5cpa = J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) - J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return 0.5 * J5cpa / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::APp_6(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J7cpa = J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false) - J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true);
    return -0.5 * J7cpa / std::sqrt(-J2ccp * J2scp);
}

double BKstarllDecay::APp_8(double q2, double m_B, double m_K, double m_l) {
    double J2ccp = J2c(q2, m_B, m_K, m_l, false) + J2c(q2, m_B, m_K, m_l, true);
    double J2scp = J2s(q2, m_B, m_K, m_l, false) + J2s(q2, m_B, m_K, m_l, true);
    double J8cpa = J8(q2, m_B, m_K, m_l, false) - J8(q2, m_B, m_K, m_l, true);
    return -J8cpa / std::sqrt(-J2ccp * J2scp);
}

std::vector<double> BKstarllDecay::dG_dq2_binned(bool bar) {
    std::vector<double> out;
    auto J_i = bar ? J_i_bar_binned : J_i_binned;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double res = 0.75 * (J_i[0][i] - (2 * J_i[1][i] + J_i[2][i]) / 3.); 
        out.push_back(res);
    }   
    return out;
}

double BKstarllDecay::dG_dq2_avg_bin(size_t bin) {
    return 0.75 * (J_i_binned[0][bin] + J_i_bar_binned[0][bin] - (2 * (J_i_binned[1][bin] + J_i_bar_binned[1][bin]) + J_i_binned[2][bin] + J_i_bar_binned[2][bin]) / 3.);
}

std::vector<double> BKstarllDecay::A_FB_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J6 = 2 * J_i_binned[7][i] + J_i_binned[9][i];
        double J6bar = 2 * J_i_bar_binned[7][i] + J_i_bar_binned[9][i];
        double res = -0.375 * (J6 + J6bar) / dG_dq2_avg_bin(i); 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_CP_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double dG = J_i_binned[0][i] - (2 * J_i_binned[1][i] + J_i_binned[2][i]) / 3.;
        double dGbar = J_i_bar_binned[0][i] - (2 * J_i_bar_binned[1][i] + J_i_bar_binned[2][i]) / 3.;
        double res = (dG - dGbar) / (dG + dGbar); 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::F_L_binned() {
    std::vector<double> out;
    std::vector<double> f_t = F_T_binned();
    for (size_t i = 0; i < this->bins.size(); i++) {
        double res = 1. - f_t[i]; 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::F_T_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double res = 4.0 * (J_i_binned[1][i] + J_i_bar_binned[1][i]) / dG_dq2_avg_bin(i); 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_1_binned(double m_B, double m_K, double m_l) {
    auto num_f = [m_B, m_K, m_l, this] (double q2) {
        double AperpApar = std::real(A_par(q2, m_B, m_K, m_l, 1, false) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, false)) + A_par(q2, m_B, m_K, m_l, -1, false) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, false)));
        double AperpApar_bar = std::real(A_par(q2, m_B, m_K, m_l, 1, true) * std::conj(A_perp(q2, m_B, m_K, m_l, 1, true)) + A_par(q2, m_B, m_K, m_l, -1, true) * std::conj(A_perp(q2, m_B, m_K, m_l, -1, true)));
        return std::pow(beta_l(q2, m_l), 2) * std::real(AperpApar + AperpApar_bar);
    };

    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scpa = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double num = integrate(num_f, this->bins[i].first, this->bins[i].second, 1e-2);
        double res = -0.5 * num / J2scpa; 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_2_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double res = 0.5 * (J_i_binned[3][i] + J_i_bar_binned[3][i]) / (J_i_binned[1][i] + J_i_bar_binned[1][i]); 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_3_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J2ccp = J_i_binned[2][i] + J_i_bar_binned[2][i];
        double J3cp = J_i_binned[3][i] + J_i_bar_binned[3][i];
        double J4cp = J_i_binned[4][i] + J_i_bar_binned[4][i];
        double J7cp = J_i_binned[11][i] + J_i_bar_binned[11][i];
        double res = std::sqrt((4 * J4cp * J4cp + J7cp * J7cp) / std::abs(2 * J2ccp * (2 * J2scp + J3cp)));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_4_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J4cp = J_i_binned[4][i] + J_i_bar_binned[4][i];
        double J5cp = J_i_binned[6][i] + J_i_bar_binned[6][i];
        double J7cp = J_i_binned[11][i] + J_i_bar_binned[11][i];
        double J8cp = J_i_binned[12][i] + J_i_bar_binned[12][i];
        double res = std::sqrt((J5cp * J5cp + 4 * J8cp * J8cp) / (J7cp * J7cp + 4 * J4cp * J4cp));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_5_binned(double m_B, double m_K, double m_l) {
    auto num_f = [m_B, m_K, m_l, this] (double q2) {
        complex_t AperpApar = A_perp(q2, m_B, m_K, m_l, 1, false) * std::conj(A_par(q2, m_B, m_K, m_l, -1, false)) + A_perp(q2, m_B, m_K, m_l, -1, false) * std::conj(A_par(q2, m_B, m_K, m_l, 1, false));
        complex_t AperpApar_bar = A_perp(q2, m_B, m_K, m_l, 1, true) * std::conj(A_par(q2, m_B, m_K, m_l, -1, true)) + A_perp(q2, m_B, m_K, m_l, -1, true) * std::conj(A_par(q2, m_B, m_K, m_l, 1, true));
        return std::pow(beta_l(q2, m_l), 2) * std::abs(AperpApar + AperpApar_bar);
    };

    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scpa = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double num = integrate(num_f, this->bins[i].first, this->bins[i].second, 1e-2);
        double res = 0.25 * num / J2scpa; 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_Re_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J6scp = J_i_binned[8][i] + J_i_bar_binned[8][i];
        double res = 0.25 * J6scp / J2scp;
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_T_Re_CPV_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J6scpa = J_i_binned[8][i] - J_i_bar_binned[8][i];
        double res = 0.25 * J6scpa / J2scp;
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::A_Im_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double res = (J_i_binned[13][i] + J_i_bar_binned[13][i]) / dG_dq2_avg_bin(i); 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::alpha_K_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J2ccp = J_i_binned[2][i] + J_i_bar_binned[2][i];
        double res = -0.5 * (2 * J2scp + J2ccp) / J2scp;
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::H_T_1_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J2ccp = J_i_binned[2][i] + J_i_bar_binned[2][i];
        double J3cp = J_i_binned[3][i] + J_i_bar_binned[3][i];
        double J4cp = J_i_binned[4][i] + J_i_bar_binned[4][i];
        double res = RT2 * J4cp / std::sqrt(std::abs(J2ccp * (J2scp - J3cp)));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::H_T_2_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J2ccp = J_i_binned[2][i] + J_i_bar_binned[2][i];
        double J3cp = J_i_binned[3][i] + J_i_bar_binned[3][i];
        double J5cp = J_i_binned[6][i] + J_i_bar_binned[6][i];
        double res = J5cp / std::sqrt(std::abs(2 * J2ccp * (2 * J2scp + J3cp)));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::H_T_3_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J3cp = J_i_binned[3][i] + J_i_bar_binned[3][i];
        double J6cp = (2 * (J_i_binned[7][i] + J_i_bar_binned[7][i]) + J_i_binned[9][i] + J_i_bar_binned[9][i]);
        double res = 0.5 * J6cp / std::sqrt(std::abs(4 * J2scp * J2scp - J3cp * J3cp));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::P_2_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J6scp = J_i_binned[7][i] + J_i_bar_binned[7][i];
        double res = 0.125 * J6scp / J2scp;
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::P_3_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J9cp = J_i_binned[13][i] + J_i_bar_binned[13][i];
        double res = -0.25 * J9cp / J2scp;
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::P_6_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J2ccp = J_i_binned[2][i] + J_i_bar_binned[2][i];
        double J3cp = J_i_binned[3][i] + J_i_bar_binned[3][i];
        double J7cp = J_i_binned[11][i] + J_i_bar_binned[11][i];
        double res = -J7cp / std::sqrt(std::abs(2 * J2ccp * (2 * J2scp - J3cp)));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::P_8_binned() {
    std::vector<double> out;
    for (size_t i = 0; i < this->bins.size(); i++) {
        double J2scp = J_i_binned[1][i] + J_i_bar_binned[1][i];
        double J2ccp = J_i_binned[2][i] + J_i_bar_binned[2][i];
        double J3cp = J_i_binned[3][i] + J_i_bar_binned[3][i];
        double J8cp = J_i_binned[12][i] + J_i_bar_binned[12][i];
        double res = -RT2 * J8cp / std::sqrt(std::abs(J2ccp * (2 * J2scp - J3cp)));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::Pp_i_binned(size_t i, bool cpv) {
    if (!(i == 4 || i == 5 || i == 6 || i == 8)) LOG_ERROR("Value Error", "P'_i(B > K*ll) is not defined for i =", i);

    std::map<size_t, double> factors = {{4, 1.0}, {5, 0.5}, {6, -0.5}, {8, -1.0}};
    std::map<size_t, size_t> J_idx = {{4, 4}, {5, 5}, {6, 10}, {8, 12}};
    double sign = cpv ? -1 : 1;

    std::vector<double> out;
    for (size_t j = 0; j < this->bins.size(); j++) {
        double J2scp = J_i_binned[1][j] + J_i_bar_binned[1][j];
        double J2ccp = J_i_binned[2][j] + J_i_bar_binned[2][j];
        double Jicp = J_i_binned[J_idx[i]][j] + sign * J_i_bar_binned[J_idx[i]][j];
        double res = Jicp / std::sqrt(std::abs(J2ccp * J2scp));
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::S_i_binned(size_t i, bool cpv) {
    if (i < 3 || i > 9) LOG_ERROR("Value Error", "S_i(B > K*ll) is not defined for i =", i);

    if (i == 6) i = 9;
    else if (i == 7) i = 10;
    else if (i >= 8) i += 4;

    double sign = cpv ? -1 : 1;

    std::vector<double> out;
    for (size_t j = 0; j < this->bins.size(); j++) {
        double res = (J_i_binned[i][j] + sign * J_i_bar_binned[i][j]) / dG_dq2_avg_bin(j); 
        out.push_back(res);
    }   
    return out;
}

std::vector<double> BKstarllDecay::P_i_CPV_binned(size_t i) {
   if (i < 1 || i > 3) LOG_ERROR("Value Error", "P_i_CPV(B > K*ll) is not defined for i =", i);

    std::map<size_t, double> factors = {{1, 0.5}, {2, 0.125}, {3, -0.25}};
    std::map<size_t, size_t> J_idx = {{1, 3}, {2, 7}, {3, 13}};

    std::vector<double> out;
    for (size_t j = 0; j < this->bins.size(); j++) {
        double J2scp = J_i_binned[1][j] + J_i_bar_binned[1][j];
        double Jicpv = J_i_binned[J_idx[i]][j] - J_i_bar_binned[J_idx[i]][j];
        double res = Jicpv / J2scp;
        out.push_back(res);
    }   
    return out;
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
    auto m_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", this->cfg.charge == Charge::B_0 ? 511 : 521));
    auto m_Ks = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", this->cfg.charge == Charge::B_0 ? 313 : 323));
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
    auto m_c_mu_b = std::make_shared<OperatorNode>("m_c(mu_b)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.m_c_mu_b = ObsQCDProxy()(MassConfig(4, w_config.hadronic_scale, MassType::POLE, MassType::POLE)); return cache.m_c_mu_b; });
    auto zc = std::make_shared<OperatorNode>("z_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.z_c = pow(cache.m_c_pole / cache.m_b_PS, 2); return cache.z_c; });
    zc->addChildren({n_m_c_pole, m_b_PS});
    // NF : For kappa I use m_b(mu_b) instead of m_b(m_b) which seems more legitimate as we use alpha_s(mu_b)
    auto kappa = std::make_shared<OperatorNode>("kappa", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { cache.kappa = 1 - 2. * cache.alpha_s_mu_b / (3. * PI) * std::log(w_config.hadronic_scale / cache.m_b_mu_b); return cache.kappa; });
    kappa->addChildren({alpha_s_mu_b, m_b_mu_b});
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
        auto bound_func = std::bind(&BKstarllDecay::T_perp_p, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8, std::placeholders::_9, std::placeholders::_10, std::placeholders::_11);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_perp_p_lookup, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], false); 
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_perp_p_bar_lookup, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], true); 
        return T_perp_p_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_perp_p_cache->addChildren({m_B, m_Ks, f_B, f_K_par, zeta_3_A, zeta_3_V, omega_10_A, delta_tilde_p, delta_tilde_m, q2_min, n_q2_high, alpha_s_mu_b, N_perp, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_perp, a_2_perp, eq, lambda_B_p, wilson_cache, wilson_bar_cache});

    auto T_perp_m_cache = std::make_shared<OperatorNode>("T_perp_m_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_perp_m, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8, std::placeholders::_9, std::placeholders::_10, std::placeholders::_11);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_perp_m_lookup, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], false);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_perp_m_bar_lookup, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], true); 
        return T_perp_m_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_perp_m_cache->addChildren({m_B, m_Ks, f_B, f_K_par, zeta_3_A, zeta_3_V, omega_10_A, delta_tilde_p, delta_tilde_m, q2_min, n_q2_high, alpha_s_mu_b, N_perp, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_perp, a_2_perp, eq, lambda_B_p, wilson_cache, wilson_bar_cache});

    auto T_par_p_cache = std::make_shared<OperatorNode>("T_par_p_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_par_p, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_par_p_lookup, values[0], values[1], false); 
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_par_p_bar_lookup, values[0], values[1], true); 
        return T_par_p_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_par_p_cache->addChildren({m_B, m_Ks, q2_min, n_q2_high, alpha_s_mu_b, N_par, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_par, a_2_par, eq, lambda_B_p, omega_0, T_par_m_0, wilson_cache, wilson_bar_cache});

    auto T_par_m_cache = std::make_shared<OperatorNode>("T_par_m_cache", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { 
        auto bound_func = std::bind(&BKstarllDecay::T_par_m, &*this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_par_m_lookup, values[0], values[1], false); 
        fill_cache(bound_func, cache.q2_min, cache.q2_high, T_par_m_bar_lookup, values[0], values[1], true); 
        return T_par_m_lookup[(size_t)(LOOKUP_SIZE / 2)]; 
    });
    T_par_m_cache->addChildren({m_B, m_Ks, q2_min, n_q2_high, alpha_s_mu_b, N_par, lambda_u_hat, m_b_PS, delta_M, zc, Lb, C_F, a_1_par, a_2_par, eq, lambda_B_p, omega_0, T_par_m_0, wilson_cache, wilson_bar_cache});

    auto binned_J_cache = std::make_shared<OperatorNode>("binned_J", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {
        compute_binned_J_i(values[0], values[1], values[2], values[4], values[5], values[3]);
        return 0; 
    });
    binned_J_cache->addChildren({m_B, m_Ks, m_l, m_s, f_B, f_K_par, q2_min, n_q2_low, n_q2_high, q2_max, m_b_PS, m_b_mu_b, m_c_mu_b, kappa, alpha_s_mu_b, C_F, N_c, lambda_B_p, N_0, tp, t0, z0, ff, T_perp_p_cache, T_perp_m_cache, T_par_m_cache});

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
            complex_t i1 = phi * (T_par_p_p_f(u, q2, m_B, m_K, false) + T_par_p_nf(u, q2, m_B, m_K, false));
            complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, m_B, m_K, false));
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
            << "," << std::real(T_perp_p_cached(q2, false)) << "," << std::imag(T_perp_p_cached(q2, false))
            << "," << std::real(T_perp_m_cached(q2, false)) << "," << std::imag(T_perp_m_cached(q2, false))
            << "," << std::real(T_par_p_cached(q2, false)) << "," << std::imag(T_par_p_cached(q2, false))
            << "," << std::real(T_par_m_cached(q2, false)) << "," << std::imag(T_par_m_cached(q2, false))
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

        LOG_INFO("q2_min =", cache.q2_min);

        std::ofstream fs;
        fs.open("B_Ksll_J.csv");
        fs << "q2,J1s,J1c,J2s,J2c,J3,J4,J5,J6s,J6c,J7,J8,J9,J1sbar,J1cbar,J2sbar,J2cbar,J3bar,J4bar,J5bar,J6sbar,J6cbar,J7bar,J8bar,J9bar\n";

        auto write_line = [&] (double q2) {
            fs << q2 
            << "," << J1s(q2, m_B, m_K, m_l, false)
            << "," << J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false)
            << "," << J2s(q2, m_B, m_K, m_l, false)
            << "," << J2c(q2, m_B, m_K, m_l, false)
            << "," << J3(q2, m_B, m_K, m_l, false)
            << "," << J4(q2, m_B, m_K, m_l, false)
            << "," << J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false)
            << "," << J6s(q2, m_B, m_K, m_l, false)
            << "," << J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false)
            << "," << J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, false)
            << "," << J8(q2, m_B, m_K, m_l, false)
            << "," << J9(q2, m_B, m_K, m_l, false)
            << "," << J1s(q2, m_B, m_K, m_l, true)
            << "," << J1c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true)
            << "," << J2s(q2, m_B, m_K, m_l, true)
            << "," << J2c(q2, m_B, m_K, m_l, true)
            << "," << J3(q2, m_B, m_K, m_l, true)
            << "," << J4(q2, m_B, m_K, m_l, true)
            << "," << J5(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true)
            << "," << J6s(q2, m_B, m_K, m_l, true)
            << "," << J6c(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true)
            << "," << J7(q2, m_B, m_K, m_l, f_B, f_K_par, m_s, true)
            << "," << J8(q2, m_B, m_K, m_l, true)
            << "," << J9(q2, m_B, m_K, m_l, true)
            << "\n";
        };

        size_t n = 200;
        double dq2 = (cache.q2_max - cache.q2_min) / n;
        double q2 = cache.q2_min;

        for (size_t i = 0; i <= n; i++) {
            write_line(q2);
            q2 += dq2;
        }

        return 0; 
    });
    test_J->addChildren({m_B, m_Ks, m_l, m_s, f_B, f_K_par, q2_min, n_q2_low, n_q2_high, q2_max, m_b_PS, m_b_mu_b, m_c_mu_b, kappa, alpha_s_mu_b, C_F, N_c, lambda_B_p, N_0, test_ff, T_perp_p_cache, T_perp_m_cache, T_par_m_cache});

    auto test_binned_obs = std::make_shared<OperatorNode>("test_binned_obs", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) {
        std::ofstream fs;
        fs.open("B_Ksll_obs.csv");
        fs << "q2_min,q2_max,dG,dGbar,afb,fl,ft,cpa,pp4,pp5,pp6,pp8\n";

        auto dG = dG_dq2_binned(false);
        auto dGbar = dG_dq2_binned(true);
        auto afb = A_FB_binned();
        auto fl = F_L_binned();
        auto ft = F_T_binned();
        auto cpa = A_CP_binned();
        auto pp4 = Pp_i_binned(4);
        auto pp5 = Pp_i_binned(5);
        auto pp6 = Pp_i_binned(6);
        auto pp8 = Pp_i_binned(8);

        auto write_line = [&] (size_t i) {
            fs << this->bins[i].first 
            << "," << this->bins[i].second 
            << "," << dG[i]
            << "," << dGbar[i]
            << "," << afb[i]
            << "," << fl[i]
            << "," << ft[i]
            << "," << cpa[i]
            << "," << pp4[i]
            << "," << pp5[i]
            << "," << pp6[i]
            << "," << pp8[i]
            << "\n";
        };

        for (size_t i = 0; i < bins.size(); i++) {
            write_line(i);
        }

        return 0; 
    });
    test_binned_obs->addChildren({binned_J_cache});

    roots.emplace(ObservableMapper::to_id(Observables::TEST_B__KS_L_L), test_binned_obs);
}