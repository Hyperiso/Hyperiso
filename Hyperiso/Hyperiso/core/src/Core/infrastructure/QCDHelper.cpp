#include "QCDHelper.h"
#include <DependentBlockManager.h>


void QCDHelper::Init() {
    LOG_DEBUG("Init QCD dependent block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double m_Z = src.at("SMINPUTS")->retrieve(4)->get_val();
        double alpha_s_mZ = src.at("SMINPUTS")->retrieve(3)->get_val();
        double m_b_mb = src.at("SMINPUTS")->retrieve(5)->get_val();
        double m_t_pole = src.at("SMINPUTS")->retrieve(6)->get_val();
        double m_c = src.at("MASS")->retrieve(4)->get_val();
        double m_s = src.at("MASS")->retrieve(3)->get_val();
        double m_d = src.at("MASS")->retrieve(2)->get_val();
        double m_u = src.at("MASS")->retrieve(1)->get_val();

        double lambda_5 = match_lambda(alpha_s_mZ, m_Z, 5);
        double lambda_6_mt_pole = match_lambda(alpha_s_explicit(m_t_pole, lambda_5, 5), m_t_pole, 6);
        double lambda_4_mb_mb = match_lambda(alpha_s_explicit(m_b_mb, lambda_5, 5), m_b_mb, 4);
        double lambda_3 = match_lambda(alpha_s_explicit(m_c, lambda_4_mb_mb, 4), m_c, 3);
        double m_t_mt = calc_mt_mt(lambda6_mt_pole, lambda_5);
        double lambda_6_mt_mt = match_lambda(alpha_s_explicit(m_t_mt, lambda_5, 5), m_t_mt, 6);
        double m_b_pole = calc_mb_pole(lambda_5);
        double lambda_4_mb_pole = match_lambda(alpha_s_explicit(m_b_pole, lambda_5, 5), m_b_pole, 4);
        double m_b_1S = calc_mb_1S(lambda_4_mb_mb, m_b_pole);
        double m_c_pole = calc_mc_pole(lambda_4_mb_mb);

        dep_block->store_or_assign(LhaID(1, 3), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(1, 3)}, lambda_3, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 4, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(1, 4, 1)}, lambda_4_mb_mb, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 4, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(1, 4, 2)}, lambda_4_mb_pole, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(1, 5)}, lambda_5, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 6, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(1, 6, 1)}, lambda_6_mt_mt, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 6, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(1, 6, 2)}, lambda_6_mt_pole, 0., 0.));
        dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", 4}, m_c_pole, 0., 0.));
        dep_block->store_or_assign(LhaID(5, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(5, 1)}, m_b_mb, 0., 0.));
        dep_block->store_or_assign(LhaID(5, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", LhaID(5, 2)}, m_b_pole, 0., 0.));
        dep_block->store_or_assign(6, std::make_shared<Parameter>(ParamId{ParameterType::SM, "QCD", 6}, m_t_mt, 0., 0.));
    };
    
    DependentBlockManager::addDependentBlock("QCD", src, ParameterType::SM, func);
}

double QCDHelper::alpha_s(double mu, MassType mass_b_type, MassType mass_t_type) {
    if (mu < (*Parameters::GetInstance())("MASS", 3)) {
        LOG_ERROR("Scale Error", "Renormalisation scale for alpha_s calculation (", mu, ") is below strange mass.");
    }
    
    return alpha_s_explicit(mu, get_lambda(mu, mass_b_type, mass_t_type), get_nf(mu, mass_b_type, mass_t_type));
}

double QCDHelper::msbar_mass(int pdg_code, double mu, MassType mass_b_type, MassType mass_t_type) {
    if (pdg_code > 6 || pdg_code < 1) {
        LOG_ERROR("ValueError", "PDG code", pdg_code, "is not a quark");
    }

    auto p = Parameters::GetInstance();
    double quark_mass = pdg_code < 5 ? (*p)("MASS", pdg_code) : pdg_code == 5 ? (*p)("QCD", LhaID(5, 1)) : (*p)("QCD", 6);
    double Qinit = pdg_code < 4 ? 1 : quark_mass;
    int n_i = get_nf(Qinit, mass_b_type, mass_t_type);
    int n_f = get_nf(mu, mass_b_type, mass_t_type);
    auto Q_bounds = getOrderedMasses(mass_b_type, mass_t_type);

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

double QCDHelper::calc_mc_pole(double lambda_4) {
    double mc = (*Parameters::GetInstance())("MASS", 4);
    double mu = (*Parameters::GetInstance())("MASS", 2);
    double md = (*Parameters::GetInstance())("MASS", 1);
    double ms = (*Parameters::GetInstance())("MASS", 3);
    double alphas_mc = alpha_s_explicit(mc, lambda_4, 4);	
 	return mc * (1 + alphas_mc / PI * (constants->C_F + alphas_mc / PI * ((13.4434 - 1.0414 * 3 
            + 1.0414 * constants->C_F * ((mu + md + ms) / mc)))));
}

double QCDHelper::calc_mb_pole(double lambda_5) {
    double mc = (*Parameters::GetInstance())("MASS", 4);
    double mu = (*Parameters::GetInstance())("MASS", 2);
    double md = (*Parameters::GetInstance())("MASS", 1);
    double ms = (*Parameters::GetInstance())("MASS", 3);
    double mb = (*Parameters::GetInstance())("SMINPUTS", 5);
    double alphas_mb = alpha_s_explicit(mb, lambda_5, 5);	
    return mb * (1. + alphas_mb / PI * (constants->C_F
            + alphas_mb / PI * ((13.4434 - 1.0414 * 4. + 1.0414 * constants->C_F * (mu + md + ms + mc) / mb))));
}

double QCDHelper::calc_mb_1S(double lambda_4, double mb_pole) {
    double mu = mb_pole / 2.;
	return mb_pole * (1 - 2. / 9 * pow(alpha_s_explicit(mu, lambda_4, 4), 2.));
}

double QCDHelper::calc_mt_mt(double lambda6_mt_pole, double lambda_5) {
    double mt_pole = (*Parameters::GetInstance())("SMINPUTS", 6);
    double alpha = alpha_s_explicit(mt_pole, lambda6_mt_pole, 6);
    double a = 307. / 32. + PI2 / 3. + PI2 / 9. * log(2.) - 1. / 6 * ZETA3 - 71. / 144. * 5.;
    double mt_mt = mt_pole / (1. + alpha / PI * (4. / 3 + alpha / PI * a));
    double lambda = match_lambda(alpha_s_explicit(mt_mt, lambda_5, 5), mt_mt, 6);
    alpha = alpha_s_explicit(mt_pole, lambda, 6);
    return mt_pole / (1. + alpha / PI * (4. / 3 + alpha / PI * a));
}

int QCDHelper::get_nf(double Q, MassType mass_b_type, MassType mass_t_type) {
    auto masses = getOrderedMasses(mass_b_type, mass_t_type);
    for (size_t i = 0; i < masses.size(); ++i) {
        if (1 - Q / masses.at(i) > 1e-4) {
            return i;
        }
    }
    return 6;
}

double QCDHelper::get_lambda(double mu, MassType mass_b_type, MassType mass_t_type) {
    int nf = get_nf(mu, mass_b_type, mass_t_type);
    if (nf == 4)
        return (*Parameters::GetInstance())("QCD", LhaID(1, 4, mass_b_type == MassType::POLE ? 2 : 1));
    if (nf == 6)
        return (*Parameters::GetInstance())("QCD", LhaID(1, 6, mass_t_type == MassType::POLE ? 2 : 1));
    return (*Parameters::GetInstance())("QCD", LhaID(1, nf));
}

std::vector<double> QCDHelper::getOrderedMasses(MassType mass_b_type, MassType mass_t_type) {
    double m_b = (*Parameters::GetInstance())("QCD", LhaID(5, mass_b_type == MassType::POLE ? 2 : 1));
    double m_t = mass_t_type == MassType::MSBAR ? (*Parameters::GetInstance())("QCD", 6) : (*Parameters::GetInstance())("SMINPUTS", 6);
    double mc = (*Parameters::GetInstance())("MASS", 4);
    double mu = (*Parameters::GetInstance())("MASS", 2);
    double md = (*Parameters::GetInstance())("MASS", 1);
    double ms = (*Parameters::GetInstance())("MASS", 3);
    return {md, mu, md, mc, m_b, m_t};
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
