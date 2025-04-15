#include "QCDHelper.h"


void QCDHelper::Init(double alpha_s_mZ, double m_Z, double mt_pole, double mb_mb, double m_c, double m_s, double m_d, double m_u) {
    param_cache.alphas_mZ = alpha_s_mZ;
    param_cache.m_Z = m_Z;
    param_cache.mb_mb = mb_mb;
    param_cache.mt_pole = mt_pole;
    param_cache.light_masses[0] = m_d;
    param_cache.light_masses[1] = m_u;
    param_cache.light_masses[2] = m_s;
    param_cache.light_masses[3] = m_c;

    lambdas_running[4] = match_lambda(param_cache.alphas_mZ, param_cache.m_Z, 5);
    lambda6_mt_pole = match_lambda(alpha_s_explicit(param_cache.mt_pole, lambdas_running[4], 5), param_cache.mt_pole, 6);
    lambdas_running[3] = match_lambda(alpha_s_explicit(param_cache.mb_mb, lambdas_running[4], 5), param_cache.mb_mb, 4);
    lambdas_running[2] = match_lambda(alpha_s_explicit(param_cache.light_masses[3], lambdas_running[3], 4), param_cache.light_masses[3], 3);

    special_masses.mt_mt = calc_mt_mt();
    lambdas_running[5] = match_lambda(alpha_s_explicit(special_masses.mt_mt, lambdas_running[4], 5), special_masses.mt_mt, 6);

    special_masses.mb_pole = calc_mb_pole();
    lambda4_mb_pole = match_lambda(alpha_s_explicit(special_masses.mb_pole, lambdas_running[4], 5), special_masses.mb_pole, 4);

    special_masses.mb_1S = calc_mb_1S();
    special_masses.mc_pole = calc_mc_pole();
}

void QCDHelper::Update() {
    auto p = Parameters::GetInstance(ParameterType::SM);
    Init(
        (*p)("SMINPUTS", 3),
        (*p)("SMINPUTS", 4),
        (*p)("SMINPUTS", 6),
        (*p)("SMINPUTS", 5),
        (*p)("MASS", 4),
        (*p)("MASS", 3),
        (*p)("MASS", 1),
        (*p)("MASS", 2)
    );
}

double QCDHelper::alpha_s(double mu, MassType mass_b_type, MassType mass_t_type) {
    if (mu < param_cache.light_masses[2]) {
        LOG_ERROR("Scale Error", "Renormalisation scale for alpha_s calculation is below strange mass(", param_cache.light_masses[2] , ").");
    }
    update_cached_values();
    set_mass_types(mass_b_type, mass_t_type);
    return alpha_s_explicit(mu, get_lambda(mu), get_nf(mu));
}

double QCDHelper::msbar_mass(int pdg_code, double mu, MassType mass_b_type, MassType mass_t_type) {
    if (pdg_code > 6 || pdg_code < 1) {
        LOG_ERROR("ValueError", "PDG code", pdg_code, "is not a quark");
    }

    update_cached_values();
    set_mass_types(mass_b_type, mass_t_type);
    double quark_mass = pdg_code < 5 ? param_cache.light_masses[pdg_code - 1] : pdg_code == 5 ? param_cache.mb_mb : special_masses.mt_mt;
    double Qinit = pdg_code < 4 ? 1 : quark_mass;
    int n_i = get_nf(Qinit);
    int n_f = get_nf(mu);
    auto Q_bounds = getOrderedMasses();

    if (fpeq(Qinit, mu))
        return quark_mass;

    while (n_i > n_f) {
        quark_mass = runMass(quark_mass, Qinit, Q_bounds.at(n_i - 1), n_i);
        Qinit = Q_bounds.at(n_i - 1);
        --n_i;
    }

    while (n_i < n_f) {
        quark_mass = runMass(quark_mass, Qinit, Q_bounds.at(n_i), n_i);
        Qinit = Q_bounds.at(n_i);
        ++n_i;
    }

    return runMass(quark_mass, Qinit, mu, n_f);
}

double QCDHelper::mass_c_pole() {
    update_cached_values();
    return special_masses.mc_pole;
}

double QCDHelper::mass_b_pole() {
    update_cached_values();
    return special_masses.mb_pole;
}

double QCDHelper::mass_b_msbar() {
    update_cached_values();
    return param_cache.mb_mb;
}

double QCDHelper::mass_b_1S() {
    update_cached_values();
    return special_masses.mb_1S;
}

double QCDHelper::mass_t_msbar() {
    update_cached_values();
    return special_masses.mt_mt;
}

double QCDHelper::mass_t_pole() {
    update_cached_values();
    return param_cache.mt_pole;
}

double QCDHelper::calc_mc_pole() {
    double alphas_mc = alpha_s_explicit(param_cache.light_masses[3], lambda4_mb_pole, 4);	
 	return param_cache.light_masses[3] * (1 + alphas_mc / PI * (constants->C_F + alphas_mc / PI * ((13.4434 - 1.0414 * 3 
            + 1.0414 * constants->C_F * ((param_cache.light_masses[0] + param_cache.light_masses[1] + param_cache.light_masses[2]) / param_cache.light_masses[3])))));
}

double QCDHelper::calc_mb_pole() {
    double alphas_mb = alpha_s_explicit(param_cache.mb_mb, lambdas_running[4], 5);	
    return param_cache.mb_mb * (1. + alphas_mb / PI * (constants->C_F
            + alphas_mb / PI * ((13.4434 - 1.0414 * 4. + 1.0414 * constants->C_F * ((param_cache.light_masses[0] + param_cache.light_masses[1] + param_cache.light_masses[2] + param_cache.light_masses[3]) / param_cache.mb_mb)))));
}

double QCDHelper::calc_mb_1S() {
    double mu = special_masses.mb_pole / 2.;
	return special_masses.mb_pole * (1 - 2. / 9 * pow(alpha_s_explicit(mu, lambda4_mb_pole, 4), 2.));
}

double QCDHelper::calc_mt_mt() {
    double alpha = alpha_s_explicit(param_cache.mt_pole, lambda6_mt_pole, 6);
    double a = 307. / 32. + PI2 / 3. + PI2 / 9. * log(2.) - 1. / 6 * ZETA3 - 71. / 144. * 5.;
    double mt_mt = param_cache.mt_pole / (1. + alpha / PI * (4. / 3 + alpha / PI * a));
    double lambda = match_lambda(alpha_s_explicit(mt_mt, lambdas_running[4], 5), mt_mt, 6);
    alpha = alpha_s_explicit(param_cache.mt_pole, lambda, 6);
    return param_cache.mt_pole / (1. + alpha / PI * (4. / 3 + alpha / PI * a));
}

void QCDHelper::update_cached_values() {
    if (!param_cache.cache_valid()) 
        QCDHelper::Update();
}

/**
 * @brief Computes the number of active flavors at a given energy scale
 * 
 * @param Q Energy scale
 * @return The number of active flavors at scale Q
 */
int QCDHelper::get_nf(double Q) {
    auto masses = getOrderedMasses();
    for (size_t i = 0; i < masses.size(); ++i) {
        if (1 - Q / masses.at(i) > 1e-4) {
            return i;
        }
    }
    return 6;
}

double QCDHelper::get_lambda(double mu) {
    int nf = get_nf(mu);
    if (nf == 4 && m_b_type == MassType::POLE)
        return lambda4_mb_pole;
    if (nf == 6 && m_t_type == MassType::POLE)
        return lambda6_mt_pole;
    return lambdas_running[nf - 1];
}

std::vector<double> QCDHelper::getOrderedMasses() {
    double m_b = m_b_type == MassType::MSBAR ? param_cache.mb_mb : special_masses.mb_pole;
    double m_t = m_t_type == MassType::MSBAR ? special_masses.mt_mt : param_cache.mt_pole;
    return {param_cache.light_masses[0], param_cache.light_masses[1], param_cache.light_masses[2], param_cache.light_masses[3], m_b, m_t};
}

double QCDHelper::match_lambda(double target_alpha, double Q, int nf) {
    auto f = [&](double L) { return alpha_s_explicit(Q, L, nf) - target_alpha; };
    double L_min = 1e-3;
    double L_max = 1.;
    double L_moy = L_min;
    double alphas_min=0;
    double alphas_max=0;
    double alphas_moy = alphas_min;

    while ((std::abs(1.- L_min / L_max) > 1e-5)&& (std::abs(1. - alphas_min / target_alpha) > 1e-4) ){
        alphas_min = alpha_s_explicit(Q, L_min, nf);
        alphas_max = alpha_s_explicit(Q, L_max, nf);
        L_moy = (L_min + L_max) / 2;
        alphas_moy = alpha_s_explicit(Q, L_moy, nf);

        (target_alpha >= alphas_min && target_alpha <= alphas_moy) ? L_max = L_moy : L_min = L_moy;
    }

    if (std::abs(1-L_min/L_max) <= 1e-5) {
        LOG_ERROR("ValueError", "Unable to find suitable QCD Lambda value to match alpha_s = " + std::to_string(target_alpha) 
                 + " at scale " + std::to_string(Q) + " GeV with " + std::to_string(nf) + " active flavors.");
        return -1;
    }

    return L_min;
}

double QCDHelper::alpha_s_explicit(double mu, double lambda, int nf) {
    double r = std::pow(mu / lambda, 2);
    double L = std::log(r);
    double LL = std::log(L);
    auto [b0, b1, b2] = constants->beta[nf - 1];
    double b02 = b0 * b0;
    double b12 = b1 * b1;
    return 4 * PI * (1 - 2 * b1 * LL / (b02 * L) + 4 * b12 * (std::pow(LL - .5, 2) + b2 * b0 / 8 / b12 - 1.25) / std::pow(b02 * L, 2)) / (b0 * L);
}

/**
 * @brief Sets the type of mass (running or pole) to be taken for the calculations
 * @param m_b_type Toggle for b mass
 * @param m_t_type Toggle for t mass
 */
void QCDHelper::set_mass_types(MassType m_b_type, MassType m_t_type) {
    if (m_b_type != QCDHelper::m_b_type)
        QCDHelper::m_b_type = m_b_type;
    if (m_t_type != QCDHelper::m_t_type)
        QCDHelper::m_t_type = m_t_type;
}

double QCDHelper::runMass(double mass, double Q_i, double Q_f, int nf) {
    return mass * R(alpha_s(Q_f, m_b_type, m_t_type), nf) 
                    / R(alpha_s(Q_i, m_b_type, m_t_type), nf);
}

double QCDHelper::R(double alpha, int nf) {
    auto [b0, b1, b2] = constants->beta[nf - 1];
    auto [g0, g1, g2] = constants->gamma[nf - 1];
    double b02 = b0 * b0;
    double a = std::pow(b0 * alpha / (2 * PI), 2 * g0 / b0);
    double b = (2 * g1 / b0 - b1 * g0 / b02) * alpha / PI;
    double c = .5 * (std::pow(2 * g1 / b0 - b1 * g0 / b02, 2) 
                   + 2 * g2 / b0 - b1 * g1 / b02 - b2 * g0 / (16 * b02) 
                   + b1 * b1 * g0 / (2 * b0 * b02)) * std::pow(alpha / PI, 2);
    return a * (1 + b + c);
}

bool QCDParamCache::cache_valid() { //TODO : bad cast
    auto p = Parameters::GetInstance(ParameterType::SM);
    bool valid = fpeq(alphas_mZ, (double)(*p)("SMINPUTS", 3)) 
                    && fpeq(m_Z, (double)(*p)("SMINPUTS", 4)) 
                    && fpeq(mb_mb, (double)(*p)("SMINPUTS", 5))
                    && fpeq(mt_pole, (double)(*p)("SMINPUTS", 6));

    for (size_t i = 0; i < 4; i++) {
        valid &= fpeq(light_masses[i], (double)(*p)("MASS", i + 1));
    }
    
    return valid;
}
