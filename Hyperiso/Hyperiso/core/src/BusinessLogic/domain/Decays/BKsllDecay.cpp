#include "BKsllDecay.h"

using Charge = BKstarllConfig::B_Charge;
using FF_Src = BKstarllConfig::FF_Src;

void BKstarllDecay::load_params() {
    fill_wilson_cache();
    fill_wilson_bar_cache();
    load_FF_params();

    ObsParameterProxy p;
    auto run = [this] (double value_1gev, double eta, double gamma) { return value_1gev * pow(eta, gamma / cache.beta_0); };
    auto gamma_perp = [this] (int n) { return 4. * cache.C_F * (psi(n + 1) + GAMMA - 1. + 1. / (n + 1)); };
    auto gamma_par = [this] (int n) { return 4. * cache.C_F * (psi(n + 2) + GAMMA - .75 - 1. /(2. * (n + 1) * (n + 2))); };

    cache.alpha_em = 1.0 / p(ParamId{ParameterType::SM, "SMINPUTS", 1});
    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_l = p(ParamId{ParameterType::SM, "MASS", 11 + 2 * (int)cfg.gen});
    cache.m_s = p(ParamId{ParameterType::SM, "MASS", 3});
    cache.mu_b = w_config.hadronic_scale;
    cache.mu_f = sqrt(cache.mu_b * p(ParamId{ParameterType::DECAY, "B_Ks", 12}));
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.alpha_s_mu_f = ObsQCDProxy()(AlphasConfig(cache.mu_f, MassType::POLE, MassType::POLE));
    cache.m_c_pole = p(ParamId{ParameterType::SM, "QCD", 4});
    cache.m_c_mu_b = ObsQCDProxy()(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_b_pole = p(ParamId{ParameterType::SM, "QCD", {5, 2}});
    cache.alpha_s_mb_pole = ObsQCDProxy()(AlphasConfig(cache.m_b_pole, MassType::POLE, MassType::POLE));
    cache.alpha_s_1_GeV = ObsQCDProxy()(AlphasConfig(1.0, MassType::POLE, MassType::POLE));
    cache.eta_f = cache.alpha_s_mu_f / cache.alpha_s_1_GeV;
    cache.m_b_mu_b = ObsQCDProxy()(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.m_b_PS = cache.m_b_pole - 4 * cache.alpha_s_mb_pole * cache.mu_f / (3 * PI);
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 511 : 521});
    cache.m_Ks = p(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 313 : 323});
    cache.zeta_3_A = p(ParamId{ParameterType::DECAY, "B_Ks", 5});
    cache.zeta_3_V = p(ParamId{ParameterType::DECAY, "B_Ks", 6});
    cache.omega_10_A = p(ParamId{ParameterType::DECAY, "B_Ks", 7});
    cache.delta_t_p = p(ParamId{ParameterType::DECAY, "B_Ks", 8});
    cache.delta_t_m = p(ParamId{ParameterType::DECAY, "B_Ks", 9});
    cache.C_F = ObsQCDProxy().get_constants()->C_F;
    cache.Nc = ObsQCDProxy().get_constants()->Nc;
    cache.beta_0 = ObsQCDProxy().get_constants()->beta[5][0]; // TODO : Link with get_nf vs. hard-coded ?
    cache.tp = std::pow(cache.m_B + cache.m_Ks, 2);
    cache.tm = std::pow(cache.m_B - cache.m_Ks, 2);
    cache.t0 = this->cfg.ff_src == BKstarllConfig::FF_Src::HLMW ? 12. : cache.tp * (1. - std::sqrt(1 - cache.tm / cache.tp));
    cache.z0 = std::real(z(0, cache.tp, cache.t0));
    cache.f_B = p(ParamId{ParameterType::FLAVOR, "FCONST", {cfg.charge == Charge::B_0 ? 511 : 521, 1}});
    cache.f_Ks_par = p(ParamId{ParameterType::FLAVOR, "FCONST", {cfg.charge == Charge::B_0 ? 313 : 323, 1}});
    cache.f_Ks_perp = run(p(ParamId{ParameterType::FLAVOR, "FCONST", {cfg.charge == Charge::B_0 ? 313 : 323, 2}}), cache.eta_f, cache.C_F);
    cache.lambda_B_p = p(ParamId{ParameterType::DECAY, "B_Ks", 10}) / (1. - cache.alpha_s_mu_f * log(pow(cache.mu_f, 2)) * 1.8 / (3. * PI));
    cache.lambda_hat_u = std::conj(p(ParamId{ParameterType::SM, "VCKM", {0, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {0, 2}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}));
    cache.a_1_perp = run(p(ParamId{ParameterType::DECAY, "B_Ks", 1}), cache.eta_f, gamma_perp(1));
    cache.a_2_perp = run(p(ParamId{ParameterType::DECAY, "B_Ks", 2}), cache.eta_f, gamma_perp(2));
    cache.a_1_par = run(p(ParamId{ParameterType::DECAY, "B_Ks", 3}), cache.eta_f, gamma_par(1));
    cache.a_2_par = run(p(ParamId{ParameterType::DECAY, "B_Ks", 4}), cache.eta_f, gamma_par(2));
    cache.omega_0 = 2. * (cache.m_B - cache.m_b_PS) / 3.;
    cache.e_q = cfg.charge == Charge::B_0 ? cache.e_d : cache.e_u;
    cache.z_c = std::pow(cache.m_c_pole / cache.m_b_PS, 2);
    cache.L_b = std::log(cache.mu_b / cache.m_b_PS);
    cache.Delta_M = -6. * cache.L_b - 4. * (1. - cache.mu_f / cache.m_b_PS);
    cache.kappa = 1 - 2. * cache.alpha_s_mu_b / (3. * PI) * std::log(cache.mu_b / cache.m_b_mu_b);
    cache.pref_T_perp = PI2 * cache.f_B * cache.f_Ks_perp / (cache.Nc * cache.m_B);
    cache.pref_T_par = PI2 * cache.f_B * cache.f_B * cache.f_Ks_par / (cache.Nc * cache.m_Ks);
    cache.N_0 = std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}) * cache.G_F * cache.alpha_em / (std::sqrt(3072. * std::pow(PI, 5) * std::pow(cache.m_B, 3)));
    cache.T_par_m_0 = 4 * cache.m_B / cache.m_b_PS * (3. * cache.lambda_hat_u * (cfg.charge == Charge::B_PLUS ? 1. : 0.) * cache.C[WCoef::C2] - cache.C_bar[WCoef::C3] - 3. * cache.C_bar[WCoef::C4]);
    cache.q2_min = 4 * std::pow(cache.m_l, 2);
    cache.q2_max = cache.tm;
    cache.q2_low = p(ParamId{ParameterType::DECAY, "B_Ksll", {10, 1}});
    cache.q2_high = p(ParamId{ParameterType::DECAY, "B_Ksll", {10, 2}});

    if (cfg.ff_type == BKstarllConfig::FF_Type::SOFT || cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) {
        auto bound_T_perp_p = std::bind(&BKstarllDecay::T_perp_p, &*this, std::placeholders::_1, std::placeholders::_2);
        fill_cache(bound_T_perp_p, cache.q2_min, cache.q2_high, cache.T_perp_p_lookup, false); 
        fill_cache(bound_T_perp_p, cache.q2_min, cache.q2_high, cache.T_perp_p_bar_lookup, true); 

        auto bound_T_perp_m = std::bind(&BKstarllDecay::T_perp_m, &*this, std::placeholders::_1, std::placeholders::_2);
        fill_cache(bound_T_perp_m, cache.q2_min, cache.q2_high, cache.T_perp_m_lookup, false); 
        fill_cache(bound_T_perp_m, cache.q2_min, cache.q2_high, cache.T_perp_m_bar_lookup, true);

        auto bound_T_par_m = std::bind(&BKstarllDecay::T_perp_m, &*this, std::placeholders::_1, std::placeholders::_2);
        fill_cache(bound_T_par_m, cache.q2_min, cache.q2_high, cache.T_par_m_lookup, false); 
        fill_cache(bound_T_par_m, cache.q2_min, cache.q2_high, cache.T_par_m_bar_lookup, true);
    }

    compute_binned_J_i();
}

void BKstarllDecay::load_FF_params() {
    int ff_id = (int)(this->cfg.ff_src) + 1;
    int sse_order = this->cfg.ff_src == FF_Src::HLMW ? 1 : 2;
    
    auto get_m = [ff_id] (int i) { return ObsParameterProxy()(ParamId{ParameterType::DECAY, "B_Ksll", {ff_id, 0, i}}); };
    cache.m_R[FF::A0] = get_m(1);
    cache.m_R[FF::V] = cache.m_R[FF::T1] = get_m(2);
    cache.m_R[FF::A1] = cache.m_R[FF::A12] = cache.m_R[FF::T2] = cache.m_R[FF::T23] = get_m(3);

    for (int i = 1; i <= 7; i++) {
        for (int j = 0; j <= sse_order; j++) {
            ParamId PId {ParameterType::DECAY, "B_Ksll", {ff_id, i, j}};
            cache.alpha_ai[(FF)(i - 1)][j] = ObsParameterProxy()(PId);
        }
    }
}

complex_t BKstarllDecay::z(double t, double t_p, double t_0) {
    double a = std::sqrt(t_p - t);
    double b = std::sqrt(t_p - t_0);
    return (a - b) / (a + b);
}

double BKstarllDecay::pole(double q2, double m_R) {
    return 1. / (1 - q2 / std::pow(m_R, 2));
}

complex_t BKstarllDecay::h(double q2, double m_q) {
    if(fpeq(m_q, 0.)) return 4./9.*(2./3.+I*PI+std::log(cache.mu_b * cache.mu_b / q2));
	
	double z=4.*m_q*m_q/q2;
    double L = 2 * std::log(m_q / cache.mu_b);
	
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

complex_t BKstarllDecay::I_1(double u, double m_q, double q2) {
    if(m_q==0.) return 1.;
	
	double epsilon=1.e-10;
    double mq2 = m_q * m_q;
    double mB2 = cache.m_B * cache.m_B;

	complex_t xp=0.5+std::sqrt(0.25-(mq2-I*epsilon)/((1.-u)*mB2+u*q2));
	complex_t xm=0.5-std::sqrt(0.25-(mq2-I*epsilon)/((1.-u)*mB2+u*q2));
	complex_t yp=0.5+std::sqrt(0.25-(mq2-I*epsilon)/q2);
	complex_t ym=0.5-std::sqrt(0.25-(mq2-I*epsilon)/q2);

	return 1.+2.*mq2/(1.-u)/(mB2-q2)*(L_1(xp)+L_1(xm)-L_1(yp)-L_1(ym));
}

complex_t BKstarllDecay::Y(double q2) {
    return h(q2, cache.m_c_pole) * (4./3. * cache.C[WCoef::C1] + cache.C[WCoef::C2] + 6. * cache.C[WCoef::C3] + 60. * cache.C[WCoef::C5])
            - 0.5 * h(q2, cache.m_b_pole) * (7. * cache.C[WCoef::C3] + 4./3. * cache.C[WCoef::C4] + 76. * cache.C[WCoef::C5] + 64./3. * cache.C[WCoef::C6])
            - 0.5 * h(q2, 0.) * (cache.C[WCoef::C3] + 4./3. * cache.C[WCoef::C4] + 16. * cache.C[WCoef::C5] + 64./3. * cache.C[WCoef::C6])
            + 4./3. * cache.C[WCoef::C3] + 64./9. * cache.C[WCoef::C5] + 64./27. * cache.C[WCoef::C6];
}

complex_t BKstarllDecay::Y_u(double q2) {
    return (h(q2, cache.m_c_pole) - h(q2, 0.)) * (4. / 3. * cache.C[WCoef::C1] + cache.C[WCoef::C2]);
}

complex_t BKstarllDecay::t_perp(double u, double m_q, double q2, double E_Kstar) {
    double mB2 = cache.m_B * cache.m_B;
    if(fpeq(q2, 0.)) {
		if (fpeq(m_q, 0.)) return 4./(1.-u);
		double epsilon=1.e-10;
		complex_t xp=0.5+std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*mB2));
		complex_t xm=0.5-std::sqrt(0.25-(m_q*m_q-I*epsilon)/((1.-u)*mB2));
		return 4./(1.-u)*(1.+2.*m_q*m_q/(1.-u)/mB2*(L_1(xp)+L_1(xm)));
	}
	else return 2.*cache.m_B/(1.-u)/E_Kstar*I_1(u,m_q,q2)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B_0((1.-u)*mB2+u*q2,m_q)-B_0(q2,m_q));
}

complex_t BKstarllDecay::t_par(double u, double m_q, double q2, double E_Kstar) {
    double mB2 = cache.m_B * cache.m_B;
    return 2.*cache.m_B/(1.-u)/E_Kstar*I_1(u,m_q,q2)+((1.-u)*mB2+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B_0((1.-u)*mB2+u*q2,m_q)-B_0(q2,m_q));
}

complex_t BKstarllDecay::F_27_u(double s_hat) {
    double z=4./s_hat;
    double l = cache.L_b;

	complex_t A=
	208./243.*l+4.*s_hat/27./(1.-s_hat)*(Li2(s_hat)+log(s_hat)*log(1.-s_hat))
	+1./729./(1.-s_hat)/(1.-s_hat)*(6.*s_hat*(29.-47.*s_hat)*log(s_hat)+785.-1600.*s_hat+833.*s_hat*s_hat+6.*PI*I*(20.-49.*s_hat+47.*s_hat*s_hat))
	-2./243./pow(1.-s_hat,3.)*(2.*std::sqrt(z-1.)*(-4.+9.*s_hat-15.*s_hat*s_hat+4.*s_hat*s_hat*s_hat)*(PI/2.-std::atan(std::sqrt(z-1.)))+9.*s_hat*s_hat*s_hat*log(s_hat)*log(s_hat)+18.*PI*I*s_hat*(1.-2.*s_hat)*log(s_hat))
	+2.*s_hat/243./pow(1.-s_hat,4.)*(36.*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)+PI2*(-4.+9.*s_hat-9.*s_hat*s_hat+3.*s_hat*s_hat*s_hat));
	
	return -6.*A;
}

complex_t BKstarllDecay::F_19_u(double s_hat) {
    double z=4./s_hat;
    double l = cache.L_b;
	
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

complex_t BKstarllDecay::F_29_u(double s_hat) {
    double z=4./s_hat;
    double l = cache.L_b;

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

complex_t BKstarllDecay::A_Seidel(double s) {
    double shat = s/std::pow(cache.m_b_PS, 2);
    double z = 4./shat;

    // In the limit s -> 0 limit
    if (fabs(s) <= 1e-6) {
    return (1.1426611796982167 + 0.517134593183505 * I) +
    shat * ((-0.3221817635475234 - 0.23271056693257725 * I) +
    shat * ((-0.16999092571433885 + 0.23271056693257727 * I) +
    ((-0.13262865040070346 + 0.6981317007977318 * I) -
    (0.11612865336448869 - 1.1635528346628865 * I) * shat) * shat)) -
    0.8559670781893004 * cache.L_b -
    shat * ((0.23868312757201646 - 0.4654211338651545 * I) +
    shat * ((-0.05761316872427984 - 0.46542113386515455 * I) +
    (-0.27983539094650206 - (0.4773662551440329 -
    0.9308422677303091 * I) * shat) * shat) +
    (-0.07407407407407407 - 0.22222222222222224 * shat) * (-shat * shat / 2));
    }

    complex_t mu_b_term = -104. / 243. * 2. * cache.L_b;
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

    return 208./243.*cache.L_b+4.*shat/27./(1.-shat)*Li2log_term
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*std::log(shat)+785.-1600.*shat+833.*shat*shat+6.*PI*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*std::sqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(PI/2.-std::atan(std::sqrt(z-1.)))+9.*shat*shat*shat*std::log(shat)*std::log(shat)+18.*PI*I*shat*(1.-2.*shat)*std::log(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*std::pow(PI/2.-std::atan(std::sqrt(z-1.)),2.)+PI2*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
}

complex_t BKstarllDecay::B_Seidel(double s) {
    double m_b = cache.m_b_PS;
    double mu_b = cache.mu_b;
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

complex_t BKstarllDecay::C_Seidel(double s) {
    double m_b = cache.m_b_PS;
    double mu_b = cache.mu_b;
    double shat = s/m_b/m_b;
    return -16./81.*std::log(s/mu_b/mu_b)+428./243.-64./27.*ZETA3+16./81.*PI*I;
}

void BKstarllDecay::fill_wilson_cache() {
    auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
    LOG_INFO(b_wilsons.size());
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

double BKstarllDecay::F_a(FF a, double q2) {
    auto ai = cache.alpha_ai.at(a);
    double P = pole(q2, cache.m_R.at(a));
    double Z = std::real(z(q2, cache.tp, cache.t0)) - cache.z0;
    return P * (ai[0] + Z * (ai[1] + Z * ai[2]));
}

double BKstarllDecay::E_K(double q2) {
    return (std::pow(cache.m_B, 2) + std::pow(cache.m_Ks, 2) - q2) / (2 * cache.m_B);
}

double BKstarllDecay::A_2(double q2) {
    double A_1 = F_a(FF::A1, q2);
    double A_12 = F_a(FF::A12, q2);
    double mB2 = cache.m_B * cache.m_B;
    double mK2 = cache.m_Ks * cache.m_Ks;
    return (cache.tp * (mB2 - mK2 - q2) * A_1 - 16. * cache.m_B * mK2 * (cache.m_B + cache.m_Ks) * A_12) / ((cache.tp - q2) * (cache.tm - q2));
}

double BKstarllDecay::T_3(double q2) {
    double T_2 = F_a(FF::T2, q2);
    double T_23 = F_a(FF::T23, q2);
    double mB2 = cache.m_B * cache.m_B;
    double mK2 = cache.m_Ks * cache.m_Ks;
    return ((mB2 - mK2) * (mB2 + 3. * mK2 - q2) * T_2 - 8. * cache.m_B * mK2 * (cache.m_B - cache.m_Ks) * T_23) / ((cache.tp - q2) * (cache.tm - q2));
}

double BKstarllDecay::xi_perp(double q2) {
    return cache.m_B * F_a(FF::V, q2) / (cache.m_B + cache.m_Ks);
}

double BKstarllDecay::xi_par(double q2) {
    return (cache.m_B + cache.m_Ks) * F_a(FF::A1, q2) / (2. * E_K(q2)) - (cache.m_B - cache.m_Ks) * A_2(q2) / cache.m_B;
}

double BKstarllDecay::f_perp(double q2) {
    return std::sqrt(2. * lambda(q2)) / (cache.m_B + cache.m_Ks) * F_a(FF::V, q2);
}

double BKstarllDecay::f_par(double q2) {
    return RT2 * (cache.m_B + cache.m_Ks) * F_a(FF::A1, q2);
}

double BKstarllDecay::f_0(double q2) {
    return ((cache.m_B * cache.m_B - q2 - cache.m_Ks * cache.m_Ks) * cache.tp * F_a(FF::A1, q2) - lambda(q2) * A_2(q2)) / (2. * cache.m_Ks * (cache.m_B + cache.m_Ks) * sqrt(q2));
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

double BKstarllDecay::gv_dga_4(double u) {
    double a1 = -60. * cache.zeta_3_A * (cache.omega_10_A + 4.) + 1680. * cache.zeta_3_V;
    double a2 = 30. * cache.zeta_3_A * (15. * cache.omega_10_A + 32.) - 12600. * cache.zeta_3_V + 36. * cache.a_1_par - 72. * cache.a_2_par - 12.;
    double a3 = -100. * cache.zeta_3_A * (9. * cache.omega_10_A + 8.) + 25200. * cache.zeta_3_V - 48. * cache.a_1_par + 240. * cache.a_2_par;
    double a4 = 525. * cache.zeta_3_A * cache.omega_10_A - 14700. * cache.zeta_3_V - 180. * cache.a_2_par;
    return -u * (a1 + u * (a2 + u * (a3 + u * a4))) / 4. + cache.delta_t_p * (9. * u - 1.5) + cache.delta_t_m * 6. * u + 3. * (cache.delta_t_p + cache.delta_t_m) * log(1 - u);
}

complex_t BKstarllDecay::F_V(double v, bool bar) {
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    return .75 * (
        h(v, cache.m_c_pole) * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6] + l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
      + h(v, cache.m_b_pole) * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] + cache.C_bar[WCoef::C6]) 
      + h(v, 0.) * (cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4] + 3. * cache.C_bar[WCoef::C6] - l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)) 
      - 8. / 27. * (cache.C_bar[WCoef::C3] - cache.C_bar[WCoef::C5] - 15. * cache.C_bar[WCoef::C6])
    );
}

double BKstarllDecay::L(double q2) {
    double mb2 = std::pow(cache.m_b_PS, 2);
    return (q2 - mb2) * std::log(1 - q2 / mb2) / q2;
}

complex_t BKstarllDecay::C_perp_0(double q2, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return C7 + q2 * (Y(q2) + cache.lambda_hat_u * Y_u(q2)) / (2. * cache.m_b_PS * cache.m_B);
}

complex_t BKstarllDecay::C_par_0(double q2, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 - cache.m_B * (Y(q2) + cache.lambda_hat_u * Y_u(q2)) / (2. * cache.m_b_PS);
}

complex_t BKstarllDecay::C_perp_f(double q2, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return C7 * (2. * std::log(cache.m_b_PS / cache.mu_b) - L(q2) + cache.Delta_M);
}

complex_t BKstarllDecay::C_par_f(double q2, double sign, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + sign * cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return -C7 * (2. * std::log(cache.m_b_PS / cache.mu_b) + 2. * L(q2) + cache.Delta_M);
}

complex_t BKstarllDecay::C_perp_nf(double q2, bool bar) {
    double s_hat = q2 / (cache.m_b_PS * cache.m_b_PS);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t F_27 = f_27(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_27_u(s_hat) * l_u;
    complex_t F_19 = f_19_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_19_u(s_hat) * l_u;
    complex_t F_29 = f_29_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_29_u(s_hat) * l_u;
    return -(
        cache.C_bar[WCoef::C2] * F_27 
      + cache.C[WCoef::C8] * f_87(s_hat, cache.L_b)
      + q2 / (2. * cache.m_b_PS * cache.m_B) * (
            (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C1] / 3.) * F_29
          + 2. * cache.C_bar[WCoef::C1] * F_19
          + cache.C[WCoef::C8] * f_89(s_hat)
        )
    ) / cache.C_F;
}

complex_t BKstarllDecay::C_par_nf(double q2, bool bar) {
    double s_hat = q2 / (cache.m_b_PS * cache.m_b_PS);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t F_27 = f_27(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_27_u(s_hat) * l_u;
    complex_t F_19 = f_19_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_19_u(s_hat) * l_u;
    complex_t F_29 = f_29_PS(s_hat, cache.L_b, cache.z_c) * (1. + l_u) + F_29_u(s_hat) * l_u;
    return (
        cache.C_bar[WCoef::C2] * F_27
      + cache.C[WCoef::C8] * f_87(s_hat, cache.L_b)
      + cache.m_B / (2 * cache.m_b_PS) * (
            (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C1] / 3.) * F_29
           + 2. * cache.C_bar[WCoef::C1] * F_19
           + cache.C[WCoef::C8] * f_89(s_hat))
    ) / cache.C_F;
}

complex_t BKstarllDecay::T_par_p_p_f(double u, double q2, bool bar) {
    return 2. * T_perp_p_p_f(u, q2, bar);
}

complex_t BKstarllDecay::T_par_p_m_f(double u, double q2, bool bar) {
    return 2. * T_perp_p_m_f(u, q2, bar);
}

complex_t BKstarllDecay::T_perp_p_p_f(double u, double q2, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] + cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return 2. * cache.m_B / (1. - u) / E_K(q2) * C7;
}

complex_t BKstarllDecay::T_perp_p_m_f(double u, double q2, bool bar) {
    complex_t C7 = cache.C[WCoef::C7] - cache.C[WCoef::CP7];
    if (bar) C7 = std::conj(C7);
    return 2. * cache.m_B / (1. - u) / E_K(q2) * C7;
} 

complex_t BKstarllDecay::T_perp_p_nf(double u, double q2, bool bar) {
    double E = E_K(q2);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t t_perp_mc = t_perp(u, cache.m_c_pole, q2, E);
    complex_t t_perp_mb = t_perp(u, cache.m_b_PS, q2, E);
    complex_t t_perp_0 = t_perp(u, 0.0, q2, E);
    return -4 * cache.e_d * cache.C[WCoef::C8] / (u + (1 - u) * q2 / (cache.m_B * cache.m_B))
            + cache.m_B / (2 * cache.m_b_PS) * (
                cache.e_u * (
                    t_perp_mc * (cache.C_bar[WCoef::C2] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] + l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.))
                    - t_perp_0 * l_u * (cache.C[WCoef::C2] - cache.C[WCoef::C1] / 6.)
                )
                + cache.e_d * (
                    t_perp_mb * (cache.C_bar[WCoef::C3] + cache.C_bar[WCoef::C4] - cache.C_bar[WCoef::C6] - 4 * cache.m_b_PS / cache.m_B * cache.C_bar[WCoef::C5])
                    + t_perp_0 * cache.C_bar[WCoef::C3]
                )
            );
}

complex_t BKstarllDecay::T_par_p_nf(double u, double q2, bool bar) {
    double E = E_K(q2);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t t_par_mc = t_par(u, cache.m_c_pole, q2, E);
    complex_t t_par_mb = t_par(u, cache.m_b_PS, q2, E);
    complex_t t_par_0 = t_par(u, 0., q2, E);
    return cache.m_B / cache.m_b_PS * (
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

complex_t BKstarllDecay::T_par_m_nf(double u, double q2, bool bar) {
    double v = cache.m_B * cache.m_B * (1 - u) + q2 * u;
    return 8. * cache.m_B * cache.m_B * cache.C[WCoef::C8] / v + 8. * cache.m_B / cache.m_b_PS * F_V(v, bar);
}

complex_t BKstarllDecay::inv_lambda_B_m(double q2) {
    double x = q2 / (cache.m_B * cache.omega_0);
    return std::exp(-x) / cache.omega_0 * (I * PI - Ei(x));
}

complex_t BKstarllDecay::I_perp_p(double q2, bool bar) {
    double pref = cache.alpha_s_mu_f / (4. * PI) * cache.C_F / cache.lambda_B_p;

    if (this->cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * (T_perp_p_p_f(u, q2, bar) + T_perp_p_nf(u, q2, bar));
        };
        return pref * c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * T_perp_p_nf(u, q2, bar);
        };
        return pref * c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BKstarllDecay::I_perp_m(double q2, bool bar) {
    double pref = cache.alpha_s_mu_f / (4. * PI) * cache.C_F / cache.lambda_B_p;

    if (this->cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * (T_perp_p_m_f(u, q2, bar) + T_perp_p_nf(u, q2, bar));
        };
        return pref * c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            return phi_Kstar(u, cache.a_1_perp, cache.a_2_perp) * T_perp_p_nf(u, q2, bar);
        };
        return pref * c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BKstarllDecay::I_par_p(double q2, bool bar) {
    if (this->cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
            double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
            complex_t i1 = phi * (T_par_p_p_f(u, q2, bar) + T_par_p_nf(u, q2, bar));
            complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
            double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
            complex_t i1 = phi * T_par_p_nf(u, q2, bar);
            complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BKstarllDecay::I_par_m(double q2, bool bar) {
    if (this->cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        auto f_soft = [q2, bar, this] (double u) {
            double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
            double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
            complex_t i1 = phi * (T_par_p_m_f(u, q2, bar) + T_par_p_nf(u, q2, bar));
            complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_soft, 0, 1, 1e-2);
    } else {
        auto f_full = [q2, bar, this] (double u) {
            double fact = cache.alpha_s_mu_f * cache.C_F / (4 * PI);
            double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
            complex_t i1 = phi * T_par_p_nf(u, q2, bar);
            complex_t i2 = phi * (cache.T_par_m_0 + fact * T_par_m_nf(u, q2, bar));
            return fact / cache.lambda_B_p * i1 + cache.e_q * inv_lambda_B_m(q2) * i2;
        };
        return c_integrate(f_full, 0, 1, 1e-2);
    }
}

complex_t BKstarllDecay::I_HSA_1(double q2, bool bar) {
    auto f = [q2, bar, this] (double u) {
        double phi = phi_Kstar(u, cache.a_1_par, cache.a_2_par);
        double v = cache.m_B * cache.m_B * (1 - u) + u * q2;
        return phi * cache.m_B * cache.m_B / v * F_V(v, bar);
    };
    
    return c_integrate(f, 0, 1, 1e-2);
}

complex_t BKstarllDecay::I_HSA_2(double q2, bool bar) {
    auto f = [q2, bar, this] (double u) {
        double int_phi_par = gv_dga_4(u);
        double v = cache.m_B * cache.m_B * (1 - u) + u * q2;
        return int_phi_par * F_V(v, bar);
    };
    
    return c_integrate(f, 0, 1, 1e-2);
}

complex_t BKstarllDecay::delta_T_perp_WA(double q2) {
    double pref = cache.e_q * 2. * PI2 * cache.f_B / (cache.m_b_PS * cache.m_B);
    complex_t W_perp = cache.C[WCoef::C3] + 4. / .3 * (cache.C[WCoef::C4] + 3. * cache.C[WCoef::C5] + 4. * cache.C[WCoef::C6]);
    complex_t W_par = cache.C_bar[WCoef::C3] + 3. * cache.C_bar[WCoef::C4];
    if (cfg.charge == Charge::B_PLUS) 
        W_par += -3. * cache.C[WCoef::C2];
    double s_hat = q2 / (cache.m_B * cache.m_B);
    return pref * (
        -2. * cache.f_Ks_perp * W_perp * F_perp(s_hat)
      + cache.f_Ks_par * cache.m_Ks * W_par / (3. * (1 - s_hat) * cache.lambda_B_p)
    );
}

complex_t BKstarllDecay::delta_T_perp_HSA(double q2, bool bar) {
    double pref = cache.e_q * cache.alpha_s_mu_b * cache.C_F * PI * cache.f_B / (cache.Nc * cache.m_b_PS * cache.m_B);
    double s_hat = q2 / (cache.m_B * cache.m_B);
    return pref * (
        3. * cache.C[WCoef::C8] * cache.m_b_PS / cache.m_B * cache.f_Ks_perp * X_perp(s_hat)
      + 2. * cache.f_Ks_perp * I_HSA_1(q2, bar)
      - cache.m_Ks * cache.f_Ks_par / ((1 - s_hat) * cache.lambda_B_p) * I_HSA_2(q2, bar)
    );
}

complex_t BKstarllDecay::T_perp_p(double q2, bool bar) {
    complex_t C_perp_p = C_perp_0(q2, 1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_perp_f(q2, 1, bar) + C_perp_nf(q2, bar));
    return xi_perp(q2) * C_perp_p + cache.pref_T_perp * I_perp_p(q2, bar) + delta_T_perp_WA(q2) + delta_T_perp_HSA(q2, bar);
}

complex_t BKstarllDecay::T_perp_m(double q2, bool bar) {
    complex_t C_perp_m = C_perp_0(q2, -1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_perp_f(q2, -1, bar) + C_perp_nf(q2, bar));
    return xi_perp(q2) * C_perp_m + cache.pref_T_perp * I_perp_m(q2, bar) + delta_T_perp_WA(q2) + delta_T_perp_HSA(q2, bar);
}

complex_t BKstarllDecay::T_par_p(double q2, bool bar) {
    complex_t C_par_p = C_par_0(q2, 1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_par_f(q2, 1, bar) + C_par_nf(q2, bar));
    return xi_par(q2) * C_par_p + cache.pref_T_par / E_K(q2) * I_par_p(q2, bar);
}

complex_t BKstarllDecay::T_par_m(double q2, bool bar) {
    complex_t C_par_m = C_par_0(q2, -1, bar) + cache.alpha_s_mu_b / (4. * PI) * (C_par_f(q2, -1, bar) + C_par_nf(q2, bar));
    return xi_par(q2) * C_par_m + cache.pref_T_par / E_K(q2) * I_par_m(q2, bar);
}

complex_t BKstarllDecay::Delta_par(double q2) {
    return 1. + cache.alpha_s_mu_b * cache.C_F / (2. * PI) * (
        L(q2) - 1 - 3. * PI2 * q2 * cache.f_B * cache.f_Ks_par * cache.m_Ks / (cache.Nc * cache.m_B * cache.lambda_B_p * xi_par(q2) * std::pow(E_K(q2), 3)) * F_perp(0.0)
    );
}

complex_t BKstarllDecay::T_perp_p_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_perp_p_bar_lookup : cache.T_perp_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_perp_m_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_perp_m_bar_lookup : cache.T_perp_m_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_par_p_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_par_p_bar_lookup : cache.T_par_p_lookup, cache.q2_min, cache.q2_high);
}

complex_t BKstarllDecay::T_par_m_cached(double q2, bool bar) {
    return lerp(q2, bar ? cache.T_par_m_bar_lookup : cache.T_par_m_lookup, cache.q2_min, cache.q2_high);
}

double BKstarllDecay::beta_l(double q2) {
    return std::sqrt(1 - 4. * cache.m_l * cache.m_l / q2);
}

double BKstarllDecay::lambda(double q2) {
    double mB2 = cache.m_B * cache.m_B;
    double mK2 = cache.m_Ks * cache.m_Ks;
    return mB2 * mB2 + mK2 * mK2 + q2 * q2 - 2. * (mB2 * mK2 + (mB2 + mK2) * q2);
}

complex_t BKstarllDecay::N(double q2, bool bar) {
    complex_t N0 = bar ? std::conj(cache.N_0) : cache.N_0; 
    return N0 * std::sqrt(q2 * beta_l(q2) * std::sqrt(lambda(q2)));
}

complex_t BKstarllDecay::delta_A_perp_QCDf(double q2, double sign, bool bar) {
    double guesstimate_err = 1.0;
    complex_t delta_A = 0.0;

    if (!fpeq(std::abs(cache.h_p_fit[0]), 0.0)) {
        complex_t h_p = cache.h_p_fit[0] + q2 * (cache.h_p_fit[1] + q2 * cache.h_p_fit[2]);
        complex_t h_m = cache.h_m_fit[0] + q2 * (cache.h_m_fit[1] + q2 * cache.h_m_fit[2]);
        delta_A = 16.0 * PI2 * RT2 * N(q2, bar) * std::pow(cache.m_B, 3) * (h_p - h_m) / q2;
    } else {
        size_t id = size_t (0.5 * (1 + sign));
        guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    }

    return 2 * RT2 * cache.m_b_PS * N(q2, bar) * std::sqrt(lambda(q2)) / q2 * T_perp_p_cached(q2, bar) * guesstimate_err + delta_A;
}

complex_t BKstarllDecay::delta_A_perp_vD(double q2, bool bar) {
    double F_perp = std::sqrt(2. * lambda(q2)) / (cache.m_B * (cache.m_B + cache.m_Ks)) * F_a(FF::V, q2);
    complex_t z_q2 = z(q2, cache.tp_nf, cache.t0_nf);
    complex_t P_H_perp = cache.alpha_perp[0] + z_q2 * (cache.alpha_perp[1] + z_q2 * cache.alpha_perp[2]);
    complex_t H_perp = (1. - z_q2 * std::conj(cache.z_Jpsi_nf)) / (z_q2 - cache.z_Jpsi_nf) * (1. - z_q2 * std::conj(cache.z_psi2S_nf)) / (z_q2 - cache.z_psi2S_nf) * P_H_perp * F_perp;
    return -32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 3) * H_perp / q2;
}

complex_t BKstarllDecay::delta_A_perp_K(double q2, bool bar) {
    double DeltaC9_M1 = (cache.r1_M[0] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[0] * cache.q2_bar / q2) / (1 + cache.r2_M[0] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    return N(q2, bar) * RT2 * std::sqrt(lambda(q2)) * DeltaC9_M1 * F_a(FF::V, q2) / (cache.m_B + cache.m_Ks);
}

complex_t BKstarllDecay::delta_A_perp(double q2, double sign, bool bar) {
    switch(cfg.power_corr_impl) {
        case BKstarllConfig::Power_Corrections_Impl::BFS:
            return delta_A_perp_QCDf(q2, sign, bar);
        case BKstarllConfig::Power_Corrections_Impl::BCvDV:
            return delta_A_perp_vD(q2, bar);
        case BKstarllConfig::Power_Corrections_Impl::KMPW:
            return delta_A_perp_K(q2, bar);
        default:
            return 0.0;
    }
}

complex_t BKstarllDecay::delta_A_par_QCDf(double q2, double sign, bool bar) {
    double guesstimate_err = 1.0;
    complex_t delta_A = 0.0;

    if (!fpeq(std::abs(cache.h_p_fit[0]), 0.0)) {
        complex_t h_p = cache.h_p_fit[0] + q2 * (cache.h_p_fit[1] + q2 * cache.h_p_fit[2]);
        complex_t h_m = cache.h_m_fit[0] + q2 * (cache.h_m_fit[1] + q2 * cache.h_m_fit[2]);
        delta_A = 16.0 * PI2 * RT2 * N(q2, bar) * std::pow(cache.m_B, 3) * (h_p + h_m) / q2;
    } else {
        size_t id = 2 + size_t (0.5 * (1 + sign));
        guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    }

    return -4 * RT2 * cache.m_b_PS * N(q2, bar) * std::sqrt(cache.tp * cache.tm) * E_K(q2) / (q2 * cache.m_B) * T_perp_m_cached(q2, bar) * guesstimate_err + delta_A;
}

complex_t BKstarllDecay::delta_A_par_vD(double q2, bool bar) {
    double F_par = RT2 * (cache.m_B + cache.m_Ks) / cache.m_B * F_a(FF::A1, q2);
    complex_t z_q2 = z(q2, cache.tp_nf, cache.t0_nf);
    complex_t P_H_par = cache.alpha_par[0] + z_q2 * (cache.alpha_par[1] + z_q2 * cache.alpha_par[2]);
    complex_t H_par = (1. - z_q2 * std::conj(cache.z_Jpsi_nf)) / (z_q2 - cache.z_Jpsi_nf) * (1. - z_q2 * std::conj(cache.z_psi2S_nf)) / (z_q2 - cache.z_psi2S_nf) * P_H_par * F_par;
    return 32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 3) * H_par / q2;
}

complex_t BKstarllDecay::delta_A_par_K(double q2, bool bar) {
    double DeltaC9_M2 = (cache.r1_M[1] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[1] * cache.q2_bar / q2) / (1 + cache.r2_M[1] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    return N(q2, bar) / RT2 * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks) * DeltaC9_M2 * F_a(FF::A1, q2) / (cache.m_B - cache.m_Ks);
}

complex_t BKstarllDecay::A_perp_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};

    if (cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * F_a(FF::V, q2) / (cache.m_B + cache.m_Ks);
        F_T = T_perp_p_cached(q2, bar);
        size_t id = size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    } else {
        F_T = (cache.C[WCoef::C7] + cache.C[WCoef::CP7]) * F_a(FF::T1, q2);
        complex_t w = cache.C[WCoef::C9] + cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] + cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) w += Y(q2);
        F = w * F_a(FF::V, q2) / (cache.m_B + cache.m_Ks);
        delta_A = delta_A_perp(q2, sign, bar);
    }

    return (N(q2, bar) * std::sqrt(2 * lambda(q2)) * (F + 2. * cache.m_b_PS * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BKstarllDecay::A_par_low(double q2, double sign, bool bar) {
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};

    if (cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * F_a(FF::A1, q2) / (cache.m_B - cache.m_Ks);
        F_T = 2. * E_K(q2) * T_perp_m_cached(q2, bar) / cache.m_B;
        size_t id = 2 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * F_a(FF::T2, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) w += Y(q2);
        F = w * F_a(FF::A1, q2) / (cache.m_B - cache.m_Ks);
        delta_A = delta_A_par(q2, sign, bar);
    }

    return (-N(q2, bar) * std::sqrt(2.) * (cache.m_B * cache.m_B - cache.m_Ks * cache.m_Ks) * (F + 2. * cache.m_b_PS * F_T / q2) + delta_A) * had_err_factor;
}

complex_t BKstarllDecay::A_0_low(double q2, double sign, bool bar) {
    double mB2 = cache.m_B * cache.m_B;
    double mK2 = cache.m_Ks * cache.m_Ks;
    complex_t F, F_T;
    complex_t delta_A {0.0};
    complex_t had_err_factor {1.0};

    if (cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) w = std::conj(w);
        F = w * 16. * cache.m_B * mK2 * F_a(FF::A12, q2);
        F_T = 2. * E_K(q2) * (mB2 + 3. * mK2 - q2) * T_perp_m_cached(q2, bar) - lambda(q2) * (T_perp_m_cached(q2, bar) + T_par_m_cached(q2, bar)) / (mB2 - mK2);
        size_t id = 4 + size_t (0.5 * (1 + sign));
        had_err_factor = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    } else {
        F_T = (cache.C[WCoef::C7] - cache.C[WCoef::CP7]) * 8. * cache.m_B * mK2 / (cache.m_B + cache.m_Ks) * F_a(FF::T23, q2);
        complex_t w = cache.C[WCoef::C9] - cache.C[WCoef::CP9] + sign * (cache.C[WCoef::C10] - cache.C[WCoef::CP10]);
        if (bar) {
            w = std::conj(w);
            F_T = std::conj(F_T);
        } 
        if (cfg.power_corr_impl == BKstarllConfig::Power_Corrections_Impl::BFS) w += Y(q2);
        F = w * 16. * cache.m_B * mK2 * F_a(FF::A12, q2);
        delta_A = delta_A_par(q2, sign, bar);
    }

    return (-N(q2, bar) / (2. * cache.m_Ks * std::sqrt(q2)) * (F + 2. * cache.m_b_PS * F_T) + delta_A) * had_err_factor;
}

complex_t BKstarllDecay::A_t_low(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }

    complex_t F;
    if (cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        F = E_K(q2) * xi_par(q2) / (cache.m_Ks * Delta_par(q2));
    } else {
        F = F_a(FF::A0, q2);
    }

    return N(q2, bar) * std::sqrt(lambda(q2) / q2) * ((C10 + q2 / (cache.m_l * (cache.m_b_mu_b + cache.m_s)) * CQ2)) * F;
}

complex_t BKstarllDecay::A_S_low(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);

    complex_t F;
    if (cfg.ff_type == BKstarllConfig::FF_Type::SOFT) {
        F = E_K(q2) * xi_par(q2) / (cache.m_Ks * Delta_par(q2));
    } else {
        F = F_a(FF::A0, q2);
    }

    return -2. * N(q2, bar) * std::sqrt(lambda(q2)) * CQ1 / (cache.m_b_mu_b + cache.m_s) * F;
}

complex_t BKstarllDecay::delta_A_par(double q2, double sign, bool bar) {
    switch(cfg.power_corr_impl) {
        case BKstarllConfig::Power_Corrections_Impl::BFS:
            return delta_A_par_QCDf(q2, sign, bar);
        case BKstarllConfig::Power_Corrections_Impl::BCvDV:
            return delta_A_par_vD(q2, bar);
        case BKstarllConfig::Power_Corrections_Impl::KMPW:
            return delta_A_par_K(q2, bar);
        default:
            return 0.0;
    }
}

complex_t BKstarllDecay::delta_A_0_QCDf(double q2, double sign, bool bar) {
    double guesstimate_err = 1.0;
    complex_t delta_A_PC = 0.0;

    if (!fpeq(std::abs(cache.h_p_fit[0]), 0.0)) {
        complex_t h_0 = cache.h_0_fit[0] + q2 * (cache.h_0_fit[1] + q2 * cache.h_0_fit[2]);
        delta_A_PC = 32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 3) * h_0 / std::sqrt(q2);
    } else {
        size_t id = 4 + size_t (0.5 * (1 + sign));
        guesstimate_err = 1.0 + cache.A_had_err_low_0[id] + cache.A_had_err_low_1[id] * q2 / 6.0;
    }

    double mB2 = cache.m_B * cache.m_B;
    double mB3 = cache.m_B * mB2;
    double mK2 = cache.m_Ks * cache.m_Ks;
    double f = lambda(q2) / ((mB2 - mK2) * mB2);
    complex_t delta_A_QCDf = -N(q2, bar) * cache.m_b_PS * mB2 / (std::sqrt(q2) * cache.m_Ks) * (
        (2 * (mB2 + 3 * mK2 - q2) * E_K(q2) / mB3 - f) * T_perp_m_cached(q2, bar)
      - f * T_par_m_cached(q2, bar)
    );
    return delta_A_QCDf * guesstimate_err + delta_A_PC;
}

complex_t BKstarllDecay::delta_A_0_vD(double q2, bool bar) {
    double F_0  = cache.m_Ks / std::sqrt(q2) * ((std::pow(cache.m_B, 2) - std::pow(cache.m_Ks, 2) - q2) * std::pow(cache.m_B + cache.m_Ks, 2) * F_a(FF::A1, q2) - lambda(q2) * A_2(q2))/(2. * std::pow(cache.m_Ks, 2) * std::pow(cache.m_B + cache.m_Ks, 2));
    complex_t z_q2 = z(q2, cache.tp_nf, cache.t0_nf);
    complex_t P_H_0 = cache.alpha_0[0] + z_q2 * cache.alpha_0[1];
    complex_t H_0 = (1. - z_q2 * std::conj(cache.z_Jpsi_nf)) / (z_q2 - cache.z_Jpsi_nf) * (1. - z_q2 * std::conj(cache.z_psi2S_nf)) / (z_q2 - cache.z_psi2S_nf) * P_H_0 * F_0;
    return 32.0 * PI2 * N(q2, bar) * std::pow(cache.m_B, 2) * (cache.m_B + cache.m_Ks) * H_0 / q2;
}

complex_t BKstarllDecay::delta_A_0_K(double q2, bool bar) {
    double DeltaC9_M2 = (cache.r1_M[1] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[1] * cache.q2_bar / q2) / (1 + cache.r2_M[1] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    double DeltaC9_M3 = (cache.r1_M[2] * (1 - cache.q2_bar / q2) + cache.DeltaC9_M_qbar[2] * cache.q2_bar / q2) / (1 + cache.r2_M[2] * (cache.q2_bar - q2) / cache.q2_Jpsi);
    return -N(q2, bar) / 2. / cache.m_Ks / std::sqrt(q2) * (((cache.m_B * cache.m_B- cache.m_Ks * cache.m_Ks - q2) * (cache.m_B + cache.m_Ks) * F_a(FF::A1, q2) * DeltaC9_M2 - lambda(q2) * A_2(q2) * DeltaC9_M3 / (cache.m_B + cache.m_Ks)));
}

complex_t BKstarllDecay::delta_A_0(double q2, double sign, bool bar) {
    switch(cfg.power_corr_impl) {
        case BKstarllConfig::Power_Corrections_Impl::BFS:
            return delta_A_0_QCDf(q2, sign, bar);
        case BKstarllConfig::Power_Corrections_Impl::BCvDV:
            return delta_A_0_vD(q2, bar);
        case BKstarllConfig::Power_Corrections_Impl::KMPW:
            return delta_A_0_K(q2, bar);
        default:
            return 0.0;
    }
}

complex_t BKstarllDecay::C7_eff(double q2, bool bar) {
    complex_t A = A_Seidel(q2);
    return (bar ? std::conj(cache.C[WCoef::C7]) : cache.C[WCoef::C7]) + cache.alpha_s_mu_b / (4. * PI) * ((cache.C[WCoef::C1]-6.*cache.C[WCoef::C2])*A-cache.C[WCoef::C8]*f_87(q2 / (cache.m_B * cache.m_B), cache.L_b));
}

complex_t BKstarllDecay::C9_eff(double q2, bool bar) {
    complex_t C_h0 = 4./3.*cache.C[WCoef::C1]+cache.C[WCoef::C2]+11./2.*cache.C[WCoef::C3]-2./3.*cache.C[WCoef::C4]+52.*cache.C[WCoef::C5]-32./3.*cache.C[WCoef::C6];
    complex_t C_hb = -0.5 * (7.*cache.C[WCoef::C3]+4./3.*cache.C[WCoef::C4]+76.*cache.C[WCoef::C5]+64./3.*cache.C[WCoef::C6]);
    complex_t C_0  = 4./3.*(cache.C[WCoef::C3]+16./3.*cache.C[WCoef::C5]+16./9.*cache.C[WCoef::C6]);
    complex_t l_u = bar ? std::conj(cache.lambda_hat_u) : cache.lambda_hat_u;
    complex_t C_mc = 8. * ((4./9.*cache.C[WCoef::C1]+1./3.*cache.C[WCoef::C2])*(1.+l_u)+2.*cache.C[WCoef::C3]+20.*cache.C[WCoef::C5]);

    complex_t A = A_Seidel(q2);
    complex_t B = B_Seidel(q2);
    complex_t C = C_Seidel(q2);

    return (bar ? std::conj(cache.C[WCoef::C9]) : cache.C[WCoef::C9])
         + h(q2, 0.) * C_h0
         + h(q2, cache.m_b_PS) * C_hb
         + C_0
         + cache.alpha_s_mu_b / (4. * PI) * (cache.C[WCoef::C1]*(B + 4. * C) - 3. * cache.C[WCoef::C2] * (2. * B - C) - cache.C[WCoef::C8] * f_89(q2 / (cache.m_B * cache.m_B)))
         + std::pow(cache.m_c_mu_b, 2) / q2 * C_mc;
}

complex_t BKstarllDecay::A_perp_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) + (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) + (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] + cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * cache.m_B / q2 * C7) * f_perp(q2) * (1 + cache.A_had_err_high[size_t (0.5 * (1 + sign))]);
}

complex_t BKstarllDecay::A_par_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * cache.m_B / q2 * C7) * f_par(q2) * (1 + cache.A_had_err_high[2 + size_t (0.5 * (1 + sign))]);
}

complex_t BKstarllDecay::A_0_high(double q2, double sign, bool bar) {
    complex_t C7 = C7_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP7]) : cache.C[WCoef::CP7]);
    complex_t C9 = C9_eff(q2, bar) - (bar ? std::conj(cache.C[WCoef::CP9]) : cache.C[WCoef::CP9]);
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    if (bar) C10 = std::conj(C10);
    return -N(q2, bar) * (C9 + sign * C10 + 2. * cache.kappa * cache.m_b_mu_b * cache.m_B / q2 * C7) * f_0(q2)  * (1 + cache.A_had_err_high[4 + size_t (0.5 * (1 + sign))]);
}

complex_t BKstarllDecay::A_t_high(double q2, bool bar) {
    complex_t C10 = cache.C[WCoef::C10] - cache.C[WCoef::CP10];
    complex_t CQ2 = cache.C[WCoef::CQ2] - cache.C[WCoef::CPQ2];
    if (bar) {
        C10 = std::conj(C10);
        CQ2 = std::conj(CQ2);
    }
    return N(q2, bar) / sqrt(q2 * lambda(q2)) * (2. * C10 + q2 / cache.m_l * CQ2 / (cache.m_b_mu_b + cache.m_s)) * F_a(FF::A0, q2) * (1 + cache.A_had_err_high[6]);
}

complex_t BKstarllDecay::A_S_high(double q2, bool bar) {
    complex_t CQ1 = cache.C[WCoef::CQ1] - cache.C[WCoef::CPQ1];
    if (bar) CQ1 = std::conj(CQ1);
    return -2. * N(q2, bar) * sqrt(lambda(q2)) * CQ1 / (cache.m_b_mu_b + cache.m_s) * F_a(FF::A0, q2) * (1 + cache.A_had_err_high[7]);
}

complex_t BKstarllDecay::interpolate(double q2, complex_t val_low, complex_t val_high) {
    if (q2 < cache.q2_low)
        return val_low;

    if (q2 > cache.q2_high)
        return val_high;

    double t = (cache.q2_high - q2) / (cache.q2_high - cache.q2_low);
    return t * val_low + (1 - t) * val_high;
}

complex_t BKstarllDecay::A_perp(double q2, double sign, bool bar) {
    return interpolate(q2, A_perp_low(q2, sign, bar), A_perp_high(q2, sign, bar));
}

complex_t BKstarllDecay::A_par(double q2, double sign, bool bar) {
    return interpolate(q2, A_par_low(q2, sign, bar), A_par_high(q2, sign, bar));
}

complex_t BKstarllDecay::A_0(double q2, double sign, bool bar) {
    return interpolate(q2, A_0_low(q2, sign, bar), A_0_high(q2, sign, bar));
}

complex_t BKstarllDecay::A_t(double q2, bool bar) {
    return interpolate(q2, A_t_low(q2, bar), A_t_high(q2, bar));
}

complex_t BKstarllDecay::A_S(double q2, bool bar) {
    return interpolate(q2, A_S_low(q2, bar), A_S_high(q2, bar));
}

double BKstarllDecay::J1s(double q2, bool bar) {
    return (2. + std::pow(beta_l(q2), 2)) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, -1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
    ) + std::pow(2. * cache.m_l, 2) / q2 * std::real(
        A_perp(q2, -1, bar) * std::conj(A_perp(q2, 1, bar))
      + A_perp(q2, -1, bar) * std::conj(A_perp(q2, 1, bar))
    );
}

double BKstarllDecay::J1c(double q2, bool bar) {
    return std::pow(std::abs(A_0(q2, -1, bar)), 2) + std::pow(std::abs(A_0(q2, 1, bar)), 2) 
         + std::pow(2 * cache.m_l, 2) / q2 * (
              std::pow(std::abs(A_t(q2, bar)), 2)
            + 2. * std::real(A_0(q2, -1, bar) * std::conj(A_0(q2, 1, bar)))
           )
         + std::pow(beta_l(q2) * std::abs(A_S(q2, bar)), 2);
}

double BKstarllDecay::J2s(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 4. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, -1, bar)), 2)
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
    );
}

double BKstarllDecay::J2c(double q2, bool bar) {
    return -std::pow(beta_l(q2), 2) * (
        std::pow(std::abs(A_0(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_0(q2, 1, bar)), 2)
    );
}

double BKstarllDecay::J3(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / 2. * (
        std::pow(std::abs(A_perp(q2, -1, bar)), 2) 
      + std::pow(std::abs(A_perp(q2, 1, bar)), 2)
      - std::pow(std::abs(A_par(q2, -1, bar)), 2)
      - std::pow(std::abs(A_par(q2, 1, bar)), 2)
    );
}

double BKstarllDecay::J4(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      + std::real(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
    );
}

double BKstarllDecay::J5(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::real(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
      - cache.m_l / std::sqrt(q2) * std::real((A_perp(q2, -1, bar) + A_perp(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BKstarllDecay::J6s(double q2, bool bar) {
    return 2. * beta_l(q2) * (
        std::real(A_perp(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::real(A_perp(q2, 1, bar) * std::conj(A_perp(q2, 1, bar))) 
    );
}

double BKstarllDecay::J6c(double q2, bool bar) {
    return 4. * beta_l(q2) * cache.m_l / std::sqrt(q2) * (std::real((A_0(q2, -1, bar) + A_0(q2, 1, bar)) * std::conj(A_S(q2, bar))));
}

double BKstarllDecay::J7(double q2, bool bar) {
    return beta_l(q2) * std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      - std::imag(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
      + cache.m_l / std::sqrt(q2) * std::imag((A_perp(q2, -1, bar) + A_perp(q2, 1, bar)) * std::conj(A_S(q2, bar)))
    );
}

double BKstarllDecay::J8(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) / std::sqrt(2.) * (
        std::imag(A_0(q2, -1, bar) * std::conj(A_perp(q2, -1, bar))) 
      + std::imag(A_0(q2, 1, bar) * std::conj(A_perp(q2, 1, bar)))
    );
}

double BKstarllDecay::J9(double q2, bool bar) {
    return std::pow(beta_l(q2), 2) * (
        std::imag(A_perp(q2, -1, bar) * std::conj(A_par(q2, -1, bar))) 
      + std::imag(A_perp(q2, 1, bar) * std::conj(A_par(q2, 1, bar))) 
    );
}

void BKstarllDecay::compute_binned_J_i() {
    auto fill_binned = [&] (std::array<std::vector<double>, 14>& dest, bool bar) {
        for (auto [q2_l, q2_u] : cfg.bins) {
            dest[0].emplace_back(integrate([&] (double q2) { return 2 * J1s(q2, bar) + J1c(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[1].emplace_back(integrate([&] (double q2) { return J2s(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[2].emplace_back(integrate([&] (double q2) { return J2c(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[3].emplace_back(integrate([&] (double q2) { return J3(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[4].emplace_back(integrate([&] (double q2) { return J4(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[5].emplace_back(integrate([&] (double q2) { return J5(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[6].emplace_back(integrate([&] (double q2) { return beta_l(q2) * J5(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[7].emplace_back(integrate([&] (double q2) { return J6s(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[8].emplace_back(integrate([&] (double q2) { return beta_l(q2) * J6s(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[9].emplace_back(integrate([&] (double q2) { return J6c(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[10].emplace_back(integrate([&] (double q2) { return J7(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[11].emplace_back(integrate([&] (double q2) { return beta_l(q2) * J7(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[12].emplace_back(integrate([&] (double q2) { return J8(q2, bar); }, q2_l, q2_u, 1e-2));
            dest[13].emplace_back(integrate([&] (double q2) { return J9(q2, bar); }, q2_l, q2_u, 1e-2));
        }
    };

    fill_binned(cache.J_i_binned, false);
    fill_binned(cache.J_i_bar_binned, true);
}

std::vector<ObservableValue> BKstarllDecay::dG_dq2_binned(bool bar) {
    std::vector<ObservableValue> out;
    auto J_i = bar ? cache.J_i_bar_binned : cache.J_i_binned;
    ObservableId id = bar ? ObservableMapper::to_id(Observables::DGAMMA_BAR_DQ2_B__KSTAR_L_L) : ObservableMapper::to_id(Observables::DGAMMA_DQ2_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.75 * (J_i[0][i] - (2 * J_i[1][i] + J_i[2][i]) / 3.); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

double BKstarllDecay::dG_dq2_avg_bin(size_t bin) {
    return 0.75 * (cache.J_i_binned[0][bin] + cache.J_i_bar_binned[0][bin] - (2 * (cache.J_i_binned[1][bin] + cache.J_i_bar_binned[1][bin]) + cache.J_i_binned[2][bin] + cache.J_i_bar_binned[2][bin]) / 3.);
}

std::vector<ObservableValue> BKstarllDecay::A_FB_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_FB_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J6 = 2 * cache.J_i_binned[7][i] + cache.J_i_binned[9][i];
        double J6bar = 2 * cache.J_i_bar_binned[7][i] + cache.J_i_bar_binned[9][i];
        double res = -0.375 * (J6 + J6bar) / dG_dq2_avg_bin(i); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

// TODO
ObservableValue BKstarllDecay::q0() {
    return ObservableValue(ObservableMapper::to_id(Observables::Q0_A_FB_B__KSTAR_L_L), 0.0);
}

std::vector<ObservableValue> BKstarllDecay::A_CP_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_CP_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double dG = cache.J_i_binned[0][i] - (2 * cache.J_i_binned[1][i] + cache.J_i_binned[2][i]) / 3.;
        double dGbar = cache.J_i_bar_binned[0][i] - (2 * cache.J_i_bar_binned[1][i] + cache.J_i_bar_binned[2][i]) / 3.;
        double res = (dG - dGbar) / (dG + dGbar); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::F_L_binned() {
    std::vector<ObservableValue> out;
    std::vector<ObservableValue> f_t = F_T_binned();
    ObservableId id = ObservableMapper::to_id(Observables::F_L_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 1. - f_t[i].value; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::F_T_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::F_T_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 4.0 * (cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i]) / dG_dq2_avg_bin(i); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_1_binned() {
    auto num_f = [this] (double q2) {
        double AperpApar = std::real(A_par(q2, 1, false) * std::conj(A_perp(q2, 1, false)) + A_par(q2, -1, false) * std::conj(A_perp(q2, -1, false)));
        double AperpApar_bar = std::real(A_par(q2, 1, true) * std::conj(A_perp(q2, 1, true)) + A_par(q2, -1, true) * std::conj(A_perp(q2, -1, true)));
        return std::pow(beta_l(q2), 2) * std::real(AperpApar + AperpApar_bar);
    };

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_1_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scpa = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double num = integrate(num_f, cfg.bins[i].first, cfg.bins[i].second, 1e-2);
        double res = -0.5 * num / J2scpa; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_2_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_2_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = 0.5 * (cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i]) / (cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i]); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_3_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_3_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J4cp = cache.J_i_binned[4][i] + cache.J_i_bar_binned[4][i];
        double J7cp = cache.J_i_binned[11][i] + cache.J_i_bar_binned[11][i];
        double res = std::sqrt((4 * J4cp * J4cp + J7cp * J7cp) / std::abs(2 * J2ccp * (2 * J2scp + J3cp)));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_4_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_4_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J4cp = cache.J_i_binned[4][i] + cache.J_i_bar_binned[4][i];
        double J5cp = cache.J_i_binned[6][i] + cache.J_i_bar_binned[6][i];
        double J7cp = cache.J_i_binned[11][i] + cache.J_i_bar_binned[11][i];
        double J8cp = cache.J_i_binned[12][i] + cache.J_i_bar_binned[12][i];
        double res = std::sqrt((J5cp * J5cp + 4 * J8cp * J8cp) / (J7cp * J7cp + 4 * J4cp * J4cp));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_5_binned() {
    auto num_f = [this] (double q2) {
        complex_t AperpApar = A_perp(q2, 1, false) * std::conj(A_par(q2, -1, false)) + A_perp(q2, -1, false) * std::conj(A_par(q2, 1, false));
        complex_t AperpApar_bar = A_perp(q2, 1, true) * std::conj(A_par(q2, -1, true)) + A_perp(q2, -1, true) * std::conj(A_par(q2, 1, true));
        return std::pow(beta_l(q2), 2) * std::abs(AperpApar + AperpApar_bar);
    };

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_5_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scpa = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double num = integrate(num_f, cfg.bins[i].first, cfg.bins[i].second, 1e-2);
        double res = 0.25 * num / J2scpa; 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_Re_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_RE_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J6scp = cache.J_i_binned[8][i] + cache.J_i_bar_binned[8][i];
        double res = 0.25 * J6scp / J2scp;
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_T_Re_CPV_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_T_RE_CPV_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J6scpa = cache.J_i_binned[8][i] - cache.J_i_bar_binned[8][i];
        double res = 0.25 * J6scpa / J2scp;
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::A_Im_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::A_IM_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double res = (cache.J_i_binned[13][i] + cache.J_i_bar_binned[13][i]) / dG_dq2_avg_bin(i); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::alpha_K_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::ALPHA_K_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double res = -0.5 * (2 * J2scp + J2ccp) / J2scp;
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::H_T_1_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::H_T_1_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J4cp = cache.J_i_binned[4][i] + cache.J_i_bar_binned[4][i];
        double res = RT2 * J4cp / std::sqrt(std::abs(J2ccp * (J2scp - J3cp)));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::H_T_2_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::H_T_2_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J5cp = cache.J_i_binned[6][i] + cache.J_i_bar_binned[6][i];
        double res = J5cp / std::sqrt(std::abs(2 * J2ccp * (2 * J2scp + J3cp)));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::H_T_3_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::H_T_3_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J6cp = (2 * (cache.J_i_binned[7][i] + cache.J_i_bar_binned[7][i]) + cache.J_i_binned[9][i] + cache.J_i_bar_binned[9][i]);
        double res = 0.5 * J6cp / std::sqrt(std::abs(4 * J2scp * J2scp - J3cp * J3cp));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_2_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_2_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J6scp = cache.J_i_binned[7][i] + cache.J_i_bar_binned[7][i];
        double res = 0.125 * J6scp / J2scp;
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_3_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_3_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J9cp = cache.J_i_binned[13][i] + cache.J_i_bar_binned[13][i];
        double res = -0.25 * J9cp / J2scp;
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_6_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_6_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J7cp = cache.J_i_binned[11][i] + cache.J_i_bar_binned[11][i];
        double res = -J7cp / std::sqrt(std::abs(2 * J2ccp * (2 * J2scp - J3cp)));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_8_binned() {
    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(Observables::P_8_B__KSTAR_L_L);
    for (size_t i = 0; i < cfg.bins.size(); i++) {
        double J2scp = cache.J_i_binned[1][i] + cache.J_i_bar_binned[1][i];
        double J2ccp = cache.J_i_binned[2][i] + cache.J_i_bar_binned[2][i];
        double J3cp = cache.J_i_binned[3][i] + cache.J_i_bar_binned[3][i];
        double J8cp = cache.J_i_binned[12][i] + cache.J_i_bar_binned[12][i];
        double res = -RT2 * J8cp / std::sqrt(std::abs(J2ccp * (2 * J2scp - J3cp)));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::Pp_i_binned(size_t i, bool cpv) {
    if (!(i == 4 || i == 5 || i == 6 || i == 8)) LOG_ERROR("Value Error", "P'_i(B > K*ll) is not defined for i =", i);

    std::map<size_t, double> factors = {{4, 1.0}, {5, 0.5}, {6, -0.5}, {8, -1.0}};
    std::map<size_t, size_t> J_idx = {{4, 4}, {5, 5}, {6, 10}, {8, 12}};
    std::map<size_t, Observables> ids = {
        {4, cpv ? Observables::P_PRIME_4_CPV_B__KSTAR_L_L : Observables::P_PRIME_4_B__KSTAR_L_L},
        {5, cpv ? Observables::P_PRIME_5_CPV_B__KSTAR_L_L : Observables::P_PRIME_5_B__KSTAR_L_L},
        {6, cpv ? Observables::P_PRIME_6_CPV_B__KSTAR_L_L : Observables::P_PRIME_6_B__KSTAR_L_L},
        {8, cpv ? Observables::P_PRIME_8_CPV_B__KSTAR_L_L : Observables::P_PRIME_8_B__KSTAR_L_L}
    };
    double sign = cpv ? -1 : 1;

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(ids[i]);
    for (size_t j = 0; j < cfg.bins.size(); j++) {
        double J2scp = cache.J_i_binned[1][j] + cache.J_i_bar_binned[1][j];
        double J2ccp = cache.J_i_binned[2][j] + cache.J_i_bar_binned[2][j];
        double Jicp = cache.J_i_binned[J_idx[i]][j] + sign * cache.J_i_bar_binned[J_idx[i]][j];
        double res = Jicp / std::sqrt(std::abs(J2ccp * J2scp));
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::S_i_binned(size_t i, bool cpv) {
    if (i < 3 || i > 9) LOG_ERROR("Value Error", "S_i(B > K*ll) is not defined for i =", i);

    if (i == 6) i = 9;
    else if (i == 7) i = 10;
    else if (i >= 8) i += 4;

    std::map<size_t, Observables> ids = {
        {3, cpv ? Observables::A_3_B__KSTAR_L_L : Observables::S_3_B__KSTAR_L_L},
        {4, cpv ? Observables::A_4_B__KSTAR_L_L : Observables::S_4_B__KSTAR_L_L},
        {5, cpv ? Observables::A_5_B__KSTAR_L_L : Observables::S_5_B__KSTAR_L_L},
        {6, cpv ? Observables::A_6S_B__KSTAR_L_L : Observables::S_6C_B__KSTAR_L_L},
        {7, cpv ? Observables::A_7_B__KSTAR_L_L : Observables::S_7_B__KSTAR_L_L},
        {8, cpv ? Observables::A_8_B__KSTAR_L_L : Observables::S_8_B__KSTAR_L_L},
        {9, cpv ? Observables::A_9_B__KSTAR_L_L : Observables::S_9_B__KSTAR_L_L}
    };

    double sign = cpv ? -1 : 1;

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(ids[i]);
    for (size_t j = 0; j < cfg.bins.size(); j++) {
        double res = (cache.J_i_binned[i][j] + sign * cache.J_i_bar_binned[i][j]) / dG_dq2_avg_bin(j); 
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

std::vector<ObservableValue> BKstarllDecay::P_i_CPV_binned(size_t i) {
   if (i < 1 || i > 3) LOG_ERROR("Value Error", "P_i_CPV(B > K*ll) is not defined for i =", i);

    std::map<size_t, double> factors = {{1, 0.5}, {2, 0.125}, {3, -0.25}};
    std::map<size_t, size_t> J_idx = {{1, 3}, {2, 7}, {3, 13}};
    std::map<size_t, Observables> ids = {
        {1, Observables::P_1_CPV_B__KSTAR_L_L},
        {2, Observables::P_2_CPV_B__KSTAR_L_L},
        {3, Observables::P_3_CPV_B__KSTAR_L_L}
    };

    std::vector<ObservableValue> out;
    ObservableId id = ObservableMapper::to_id(ids[i]);
    for (size_t j = 0; j < cfg.bins.size(); j++) {
        double J2scp = cache.J_i_binned[1][j] + cache.J_i_bar_binned[1][j];
        double Jicpv = cache.J_i_binned[J_idx[i]][j] - cache.J_i_bar_binned[J_idx[i]][j];
        double res = Jicpv / J2scp;
        out.emplace_back(id, res, cfg.bins[i]);
    }   
    return out;
}

void BKstarllDecay::test_ff() {
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
}

void BKstarllDecay::test_T() {
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
    double dq2 = (cache.q2_high - cache.q2_min) / n;
    double q2 = cache.q2_min;
    for (size_t i = 0; i <= n; i++) {
        write_line(q2);
        q2 += dq2;
    }
}

void BKstarllDecay::test_J() {
    std::ofstream fs;
    fs.open("B_Ksll_J.csv");
    fs << "q2,J1s,J1c,J2s,J2c,J3,J4,J5,J6s,J6c,J7,J8,J9,J1sbar,J1cbar,J2sbar,J2cbar,J3bar,J4bar,J5bar,J6sbar,J6cbar,J7bar,J8bar,J9bar\n";

    auto write_line = [&] (double q2) {
        fs << q2 
        << "," << J1s(q2, false)
        << "," << J1c(q2, false)
        << "," << J2s(q2, false)
        << "," << J2c(q2, false)
        << "," << J3(q2, false)
        << "," << J4(q2, false)
        << "," << J5(q2, false)
        << "," << J6s(q2, false)
        << "," << J6c(q2, false)
        << "," << J7(q2, false)
        << "," << J8(q2, false)
        << "," << J9(q2, false)
        << "," << J1s(q2, true)
        << "," << J1c(q2, true)
        << "," << J2s(q2, true)
        << "," << J2c(q2, true)
        << "," << J3(q2, true)
        << "," << J4(q2, true)
        << "," << J5(q2, true)
        << "," << J6s(q2, true)
        << "," << J6c(q2, true)
        << "," << J7(q2, true)
        << "," << J8(q2, true)
        << "," << J9(q2, true)
        << "\n";
    };

    size_t n = 200;
    double dq2 = (cache.q2_max - cache.q2_min) / n;
    double q2 = cache.q2_min;

    for (size_t i = 0; i <= n; i++) {
        write_line(q2);
        q2 += dq2;
    }
}

void BKstarllDecay::test_binned_obs() {
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
        fs << cfg.bins[i].first 
        << "," << cfg.bins[i].second 
        << "," << dG[i].value
        << "," << dGbar[i].value
        << "," << afb[i].value
        << "," << fl[i].value
        << "," << ft[i].value
        << "," << cpa[i].value
        << "," << pp4[i].value
        << "," << pp5[i].value
        << "," << pp6[i].value
        << "," << pp8[i].value
        << "\n";
    };

    for (size_t i = 0; i < cfg.bins.size(); i++) {
        write_line(i);
    }
}

std::vector<ObservableValue> BKstarllDecay::compute_observable(Observables obs) {
    switch (obs) {
    case Observables::TEST:   
        test_binned_obs();
        return {};
    case Observables::DGAMMA_DQ2_B__KSTAR_L_L:   
        return dG_dq2_binned(false);
    case Observables::DGAMMA_BAR_DQ2_B__KSTAR_L_L:   
        return dG_dq2_binned(true);
    case Observables::A_FB_B__KSTAR_L_L:   
        return A_FB_binned();
    case Observables::Q0_A_FB_B__KSTAR_L_L:   
        return {q0()};
    case Observables::A_CP_B__KSTAR_L_L:   
        return A_CP_binned();
    case Observables::F_L_B__KSTAR_L_L:   
        return F_L_binned();
    case Observables::F_T_B__KSTAR_L_L:   
        return F_T_binned();
    case Observables::A_T_1_B__KSTAR_L_L:   
        return A_T_1_binned();
    case Observables::A_T_2_B__KSTAR_L_L:   
        return A_T_2_binned();
    case Observables::A_T_3_B__KSTAR_L_L:   
        return A_T_3_binned();
    case Observables::A_T_4_B__KSTAR_L_L:   
        return A_T_4_binned();
    case Observables::A_T_5_B__KSTAR_L_L:   
        return A_T_5_binned();
    case Observables::A_T_RE_B__KSTAR_L_L:   
        return A_T_Re_binned();
    case Observables::A_T_RE_CPV_B__KSTAR_L_L:   
        return A_T_Re_CPV_binned();
    case Observables::A_IM_B__KSTAR_L_L:   
        return A_Im_binned();
    case Observables::ALPHA_K_B__KSTAR_L_L:   
        return alpha_K_binned();
    case Observables::H_T_1_B__KSTAR_L_L:   
        return H_T_1_binned();
    case Observables::H_T_2_B__KSTAR_L_L:   
        return H_T_2_binned();
    case Observables::H_T_3_B__KSTAR_L_L:   
        return H_T_3_binned();
    case Observables::P_2_B__KSTAR_L_L:   
        return P_2_binned();
    case Observables::P_3_B__KSTAR_L_L:   
        return P_3_binned();
    case Observables::P_6_B__KSTAR_L_L:   
        return P_6_binned();
    case Observables::P_8_B__KSTAR_L_L:   
        return P_8_binned();
    case Observables::P_PRIME_4_B__KSTAR_L_L:   
        return Pp_i_binned(4);
    case Observables::P_PRIME_5_B__KSTAR_L_L:   
        return Pp_i_binned(5);
    case Observables::P_PRIME_6_B__KSTAR_L_L:   
        return Pp_i_binned(6);
    case Observables::P_PRIME_8_B__KSTAR_L_L:   
        return Pp_i_binned(8);
    case Observables::S_3_B__KSTAR_L_L:   
        return S_i_binned(3);
    case Observables::S_4_B__KSTAR_L_L:   
        return S_i_binned(4);
    case Observables::S_5_B__KSTAR_L_L:   
        return S_i_binned(5);
    case Observables::S_6C_B__KSTAR_L_L:   
        return S_i_binned(6);
    case Observables::S_7_B__KSTAR_L_L:   
        return S_i_binned(7);
    case Observables::S_8_B__KSTAR_L_L:   
        return S_i_binned(8);
    case Observables::S_9_B__KSTAR_L_L:   
        return S_i_binned(9);
    case Observables::A_3_B__KSTAR_L_L:   
        return S_i_binned(3, true);
    case Observables::A_4_B__KSTAR_L_L:   
        return S_i_binned(4, true);
    case Observables::A_5_B__KSTAR_L_L:   
        return S_i_binned(5, true);
    case Observables::A_6S_B__KSTAR_L_L:   
        return S_i_binned(6, true);
    case Observables::A_7_B__KSTAR_L_L:   
        return S_i_binned(7, true);
    case Observables::A_8_B__KSTAR_L_L:   
        return S_i_binned(8, true);
    case Observables::A_9_B__KSTAR_L_L:   
        return S_i_binned(9, true);
    case Observables::P_1_CPV_B__KSTAR_L_L:   
        return P_i_CPV_binned(1);
    case Observables::P_2_CPV_B__KSTAR_L_L:   
        return P_i_CPV_binned(2);
    case Observables::P_3_CPV_B__KSTAR_L_L:   
        return P_i_CPV_binned(3);
    case Observables::P_PRIME_4_CPV_B__KSTAR_L_L:   
        return Pp_i_binned(4, true);
    case Observables::P_PRIME_5_CPV_B__KSTAR_L_L:   
        return Pp_i_binned(5, true);
    case Observables::P_PRIME_6_CPV_B__KSTAR_L_L:   
        return Pp_i_binned(6, true);
    case Observables::P_PRIME_8_CPV_B__KSTAR_L_L:   
        return Pp_i_binned(8, true);
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }
}

std::vector<ObservableValue> BKstarllDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
