#include "EWHelper.h"

void EWHelper::Init() {
    LOG_DEBUG("Init EW dependent block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::SM, {"SMINPUTS", "QCD"}},
        {ParameterType::WILSON, {"B_SCALE", "D_SCALE"}}
    };

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        double inv_alpha_m_Z = src.get_val("SMINPUTS", 1);
        double m_Z = src.get_val("SMINPUTS", 4);
        double m_b_pole = src.get_val("QCD", {5, 2});
        double mu_b = src.get_val("B_SCALE", 1);
        double mu_c = src.get_val("D_SCALE", 1);

        double alpha_Z = 1 / inv_alpha_m_Z;
        double alpha_b = alpha_5_mu_b(inv_alpha_m_Z, m_Z, mu_b);
        double alpha_c = alpha_4_mu_c(inv_alpha_m_Z, m_Z, m_b_pole, mu_b, mu_c);
        double alpha_0 = 1 / 137.036; // MAJ : Store somewhere

        dep_block->store_or_assign({1, 1}, std::make_shared<Parameter>(ParamId{ParameterType::SM, "EW", {1, 1}}, alpha_Z, 0., 0.));
        dep_block->store_or_assign({1, 2}, std::make_shared<Parameter>(ParamId{ParameterType::SM, "EW", {1, 2}}, alpha_b, 0., 0.));
        dep_block->store_or_assign({1, 3}, std::make_shared<Parameter>(ParamId{ParameterType::SM, "EW", {1, 3}}, alpha_c, 0., 0.));
        dep_block->store_or_assign({1, 4}, std::make_shared<Parameter>(ParamId{ParameterType::SM, "EW", {1, 4}}, alpha_0, 0., 0.));
    };
    
    DependentBlockManager::addDependentBlock("EW", src, ParameterType::SM, func);
}

double EWHelper::alpha_em(double) {
    LOG_ERROR("NotImplementedError", "Running of alpha_em is not implemented yet. Use special values in EW block.");
    return 0.0;
}

double EWHelper::I_s(int k, double mu_low, double mu_high) {
    auto f = [k] (double s) {
        return QCDHelper::alpha_s(std::sqrt(s)) / s;
    }; 

    return integrate(f, std::pow(mu_low, 2), std::pow(mu_high, 2), 1e-3);
}

double EWHelper::delta_lept(double mu_low, double mu_high, double alpha) {
    return 3 * (1 + 3. / 4 * alpha / PI) * std::log(std::pow(mu_high / mu_low, 2));
}

double EWHelper::delta_part(double mu_low, double mu_high, int n_f, double alpha) {
    double a = 2. + n_f /3.0;
    double b = 5. / 6 + n_f / 36.0;
    return (a + b * alpha / PI) * std::log(std::pow(mu_high / mu_low, 2));
}

double EWHelper::delta_had(double mu_low, double mu_high, int n_f, double alpha) {
    double a_1 = n_f == 4 ? 10. / 3 - 17. / 34 * alpha / PI : 11. / 3 - 35. / 108 * alpha / PI;
    double a_2 = n_f == 4 ? 10. / 3 * 287. / 144 : 11. / 3 * 265. / 144;
    double a_3 = n_f == 4 ? 200675. / 23328 - 335. / 81 * ZETA3 : 257543. / 46656 - 620. / 81 * ZETA3;

    return a_1 * I_s(1, mu_low, mu_high) + a_2 * I_s(2, mu_low, mu_high) + a_3 * I_s(3, mu_low, mu_high);
}

double EWHelper::pi_heavy(double mu, double e_q, double m_q_pole, double n_l, double alpha) {
    double l_mu = std::log(std::pow(mu / m_q_pole, 2));
    double a_s = QCDHelper::alpha_s(mu) / PI;
    double pi_0_0 = l_mu;
    double pi_1_0 = e_q * e_q * (45. / 16 + 3. / 4 * l_mu);
    double pi_0_1 = 15. / 4 + l_mu;
    double pi_0_2 = 41219. / 2592 - 917. / 1296 * n_l + (4 + 4. / 3 * std::log(2) - 2. / 3 * n_l) * PI2 / 6 + 607. / 144 * ZETA3 + (437. / 36 - 7. / 9 * n_l) * l_mu + (31. / 24 - 1. / 12 * n_l) * l_mu * l_mu;

    return pi_0_0 + alpha / PI * pi_1_0 + a_s * pi_0_1 + a_s * a_s * pi_0_2;
}

double EWHelper::pi_light_heavy(double mu, double m_q_pole, double n_l) {
    double l_mu = std::log(std::pow(mu / m_q_pole, 2));
    double a_s = QCDHelper::alpha_s(mu) / PI;
    double Q_factor = 2 * (4. / 9) + (n_l - 2) * (1. / 9);
    return a_s * a_s * 3 * Q_factor * (295. / 1296 - 11. / 72 * l_mu - 1. / 12 * l_mu * l_mu);
}

double EWHelper::alpha_5_mu_b(double inv_alpha_m_Z, double m_Z, double mu_b) {
    double delta_lept_bZ = delta_lept(mu_b, m_Z, 1 / inv_alpha_m_Z);
    double delta_part_bZ = delta_part(mu_b, m_Z, 5, 1 / inv_alpha_m_Z);
    double delta_had_bZ = delta_had(mu_b, m_Z, 5, 1 / inv_alpha_m_Z);
    double inv_alpha_b = inv_alpha_m_Z + 1. / (3 * PI) * (delta_lept_bZ + delta_part_bZ + delta_had_bZ);
    return 1 / inv_alpha_b;
}

double EWHelper::alpha_4_mu_b(double inv_alpha_m_Z, double m_Z, double m_b_pole, double mu_b) {
    double alpha_5 = alpha_5_mu_b(inv_alpha_m_Z, m_Z, mu_b);
    double delta_b = 4 / 3 * pi_heavy(mu_b, -2. / 3, m_b_pole, 4, 1 / inv_alpha_m_Z) + pi_light_heavy(mu_b, m_b_pole, 4);
    double inv_alpha_4 = 1 / alpha_5 + delta_b / (3 * PI);
    return 1 / inv_alpha_4;
}

double EWHelper::alpha_4_mu_c(double inv_alpha_m_Z, double m_Z, double m_b_pole, double mu_b, double mu_c) {
    double alpha_b = alpha_4_mu_b(inv_alpha_m_Z, m_Z, m_b_pole, mu_b);
    double delta_lept_cb = delta_lept(mu_c, mu_b, 1 / inv_alpha_m_Z);
    double delta_part_cb = delta_part(mu_c, mu_b, 4, 1 / inv_alpha_m_Z);
    double delta_had_cb = delta_had(mu_c, mu_b, 4, 1 / inv_alpha_m_Z);
    double inv_alpha_mu_c = 1 / alpha_b + 1. / (3 * PI) * (delta_lept_cb + delta_part_cb + delta_had_cb);
    return 1 / inv_alpha_mu_c;
}
