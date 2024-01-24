#include "epsilon_calculator.h"
#include "../Math/Math.h"

#include <cmath>

// Constructeur de Parameters
Parameters::Parameters() {
    // Initialisation des paramètres utilisés dans les calculs epsilon_x
    SM = 1;  // Exemple de valeur, à ajuster selon votre cas
    gp = 0;
    g2 = 0;
    MSOFT_Q = 0;
    mass_top_pole = 0;
    mass_b_pole = 0;
    A_b = 0;
    tan_beta = 1;  // Éviter la division par zéro, ajustez selon votre cas
    mu_Q = 0;
    mass_gluino = 1; // Éviter la division par zéro
    mass_b1 = 0;
    mass_b2 = 0;
    inv_alpha_em = 137; // Valeur exemple pour la constante de structure fine
    M2_Q = 0;
    sbot_mix = std::vector<std::vector<double>>(3, std::vector<double>(3, 0)); // sbot_mix[1][1], sbot_mix[1][2], etc.
    mass_t1 = 0;
    mass_t2 = 0;
    yut = std::vector<double>(4, 0); // yut[3]
    charg_Umix = std::vector<std::vector<double>>(3, std::vector<double>(3, 0)); // charg_Umix[1][2], charg_Umix[2][2], etc.
    charg_Vmix = std::vector<std::vector<double>>(3, std::vector<double>(3, 0)); // charg_Vmix[1][2], charg_Vmix[2][2], etc.
    mass_cha1 = 1; // Éviter la division par zéro
    mass_cha2 = 1; // Éviter la division par zéro
    stop_mix = std::vector<std::vector<double>>(3, std::vector<double>(3, 0)); // stop_mix[1][1], stop_mix[1][2], etc.
    neut_mix = std::vector<std::vector<double>>(6, std::vector<double>(5, 0)); // neut_mix[ie][4], neut_mix[ie][3], etc.
    mass_neut = std::vector<double>(6, 0); // mass_neut[5]
    yub = std::vector<double>(4, 0); // yub[3]
    stop_tan_betamix = std::vector<std::vector<double>>(2, std::vector<double>(2, 0.0));

    A_t = 0;
    MqL3_Q = 0;
    MbR_Q = 0;
    mass_stl = 0;
}


EpsilonCalculator::EpsilonCalculator(const Parameters& p) : param(p){}


double EpsilonCalculator::epsilon_0() {
    if(param.SM == 1) return 0.;

    double sw2 = std::pow(std::sin(std::atan(param.gp / param.g2)), 2);
    double alphas_MSOFT = param.run.runningAlphasCalculation(param.MSOFT_Q);

    // Supposons que les valeurs suivantes sont définies dans param
    // Exemple: param.A_b, param.tan_beta, param.mu_Q, param.mass_gluino, etc.

    double term1 = 2.0 / 3.0 * alphas_MSOFT / M_PI * ((param.A_b / param.tan_beta - param.mu_Q) / param.mass_gluino *
               H2(param.mass_b1 * param.mass_b1 / param.mass_gluino / param.mass_gluino, param.mass_b2 * param.mass_b2 / param.mass_gluino / param.mass_gluino));
    double term2 = -0.5 * (B(param.mass_gluino, param.mass_b1, param.MSOFT_Q) + B(param.mass_gluino, param.mass_b2, param.MSOFT_Q)) / param.tan_beta;
    double term3 = 1.0 / param.inv_alpha_em / sw2 / 4.0 / M_PI * (param.mu_Q * param.M2_Q) * 
               ((param.sbot_mix[1][1] * param.sbot_mix[1][1] * H2(param.M2_Q * param.M2_Q / param.mass_b1 / param.mass_b1, param.mu_Q * param.mu_Q / param.mass_b1 / param.mass_b1) / param.mass_b1 / param.mass_b1 / 2.0) +
               (param.sbot_mix[1][2] * param.sbot_mix[1][2] * H2(param.M2_Q * param.M2_Q / param.mass_b2 / param.mass_b2, param.mu_Q * param.mu_Q / param.mass_b2 / param.mass_b2) / param.mass_b2 / param.mass_b2 / 2.0));

    return term1 + term2 + term3;
}


double EpsilonCalculator::epsilon_2() const {
    if(param.SM == 1) return 0.;

    double sw2 = std::pow(std::sin(std::atan(param.gp / param.g2)), 2);

    // Supposons que les valeurs suivantes sont définies dans param
    // Exemple: param.yut[3], param.mu_Q, param.tan_beta, param.A_t, etc.

    double term1 = param.yut[3] * param.yut[3] / 16.0 / M_PI / M_PI * 
                   (param.mu_Q / param.tan_beta - param.A_t) * 
                   ((param.charg_Umix[1][2] * param.charg_Vmix[1][2] / param.mass_cha1 * 
                     H2(param.mass_t1 * param.mass_t1 / param.mass_cha1 / param.mass_cha1, param.mass_t2 * param.mass_t2 / param.mass_cha1 / param.mass_cha1)) +
                    (param.charg_Umix[2][2] * param.charg_Vmix[2][2] / param.mass_cha2 * 
                     H2(param.mass_t1 * param.mass_t1 / param.mass_cha2 / param.mass_cha2, param.mass_t2 * param.mass_t2 / param.mass_cha2 / param.mass_cha2)));
    
    double term2 = 1.0 / param.inv_alpha_em / sw2 / 4.0 / M_PI * (param.mu_Q * param.M2_Q) * 
                   ((param.stop_mix[1][1] * param.stop_mix[1][1] * 
                     H2(param.M2_Q * param.M2_Q / param.mass_t1 / param.mass_t1, param.mu_Q * param.mu_Q / param.mass_t1 / param.mass_t1) / param.mass_t1 / param.mass_t1) +
                    (param.stop_mix[1][2] * param.stop_mix[1][2] * 
                     H2(param.M2_Q * param.M2_Q / param.mass_t2 / param.mass_t2, param.mu_Q * param.mu_Q / param.mass_t2 / param.mass_t2) / param.mass_t2 / param.mass_t2));

    return term1 + term2;
}

double EpsilonCalculator::epsilon_b() {
    return epsilon_0() + epsilon_2();
}

// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_bp
double EpsilonCalculator::epsilon_bp() {
    if(param.SM == 1) return 0.;

    double sw2 = std::pow(std::sin(std::atan(param.gp / param.g2)), 2);
    double alphas_MSOFT = param.run.runningAlphasCalculation(param.MSOFT_Q);
    int nb_neut = (param.mass_neut[5] == 0.) ? 4 : 5;

    double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       (param.A_b / param.tan_beta - param.mu_Q) / param.mass_gluino *
                       (param.stop_mix[1][1] * param.stop_mix[1][1] * param.sbot_mix[1][1] * param.sbot_mix[1][1] *
                        H2(param.mass_t1 * param.mass_t1 / param.mass_gluino / param.mass_gluino, param.mass_b2 * param.mass_b2 / param.mass_gluino / param.mass_gluino) +
                        param.stop_mix[1][1] * param.stop_mix[1][1] * param.sbot_mix[1][2] * param.sbot_mix[1][2] *
                        H2(param.mass_t1 * param.mass_t1 / param.mass_gluino / param.mass_gluino, param.mass_b1 * param.mass_b1 / param.mass_gluino / param.mass_gluino) +
                        param.stop_mix[1][2] * param.stop_tan_betamix[1][2] * param.sbot_mix[1][1] * param.sbot_mix[1][1] *
                        H2(param.mass_t2 * param.mass_t2 / param.mass_gluino / param.mass_gluino, param.mass_b2 * param.mass_b2 / param.mass_gluino / param.mass_gluino) +
                        param.stop_mix[1][2] * param.stop_mix[1][2] * param.sbot_mix[1][2] * param.sbot_mix[1][2] *
                        H2(param.mass_t2 * param.mass_t2 / param.mass_gluino / param.mass_gluino, param.mass_b1 * param.mass_b1 / param.mass_gluino / param.mass_gluino));

    for(int ie = 1; ie <= nb_neut; ++ie) {
        epsilonbp += param.yut[3] * param.yut[3] / 16.0 / M_PI / M_PI * 
                     param.neut_mix[ie][4] * param.neut_mix[ie][3] * 
                     (param.A_t - param.mu_Q / param.tan_beta) / param.mass_neut[ie] *
                     (param.stop_mix[1][1] * param.stop_mix[1][1] * param.sbot_mix[1][1] * param.sbot_mix[1][1] *
                      H2(param.mass_t2 * param.mass_t2 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b1 * param.mass_b1 / param.mass_neut[ie] / param.mass_neut[ie]) +
                      param.stop_mix[1][1] * param.stop_mix[1][1] * param.sbot_mix[1][2] * param.sbot_mix[1][2] *
                      H2(param.mass_t2 * param.mass_t2 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b2 * param.mass_b2 / param.mass_neut[ie] / param.mass_neut[ie]) +
                      param.stop_mix[1][2] * param.stop_mix[1][2] * param.sbot_mix[1][1] * param.sbot_mix[1][1] *
                      H2(param.mass_t1 * param.mass_t1 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b1 * param.mass_b1 / param.mass_neut[ie] / param.mass_neut[ie]) +
                      param.stop_mix[1][2] * param.stop_mix[1][2] * param.sbot_mix[1][2] * param.sbot_mix[1][2] *
                      H2(param.mass_t1 * param.mass_t1 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b2 * param.mass_b2 / param.mass_neut[ie] / param.mass_neut[ie]));
    }

    epsilonbp += 1.0 / param.inv_alpha_em / sw2 / 4.0 / M_PI * 
                 (param.mu_Q * param.M2_Q) * 
                 ((param.stop_mix[1][1] * param.stop_mix[1][1] *
                   H2(param.M2_Q * param.M2_Q / param.mass_t1 / param.mass_t1, param.mu_Q * param.mu_Q / param.mass_t1 / param.mass_t1) / param.mass_t1 / param.mass_t1 +
                   param.stop_mix[1][2] * param.stop_mix[1][2] *
                   H2(param.M2_Q * param.M2_Q / param.mass_t2 / param.mass_t2, param.mu_Q * param.mu_Q / param.mass_t2 / param.mass_t2) / param.mass_t2 / param.mass_t2) / 2.0 +
                 (param.sbot_mix[1][1] * param.sbot_mix[1][1] *
                  H2(param.M2_Q * param.M2_Q / param.mass_b1 / param.mass_b1, param.mu_Q * param.mu_Q / param.mass_b1 / param.mass_b1) / param.mass_b1 / param.mass_b1 +
                  param.sbot_mix[1][2] * param.sbot_mix[1][2] *
                  H2(param.M2_Q * param.M2_Q / param.mass_b2 / param.mass_b2, param.mu_Q * param.mu_Q / param.mass_b2 / param.mass_b2) / param.mass_b2 / param.mass_b2));

    return epsilonbp;
}

// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_0p
double EpsilonCalculator::epsilon_0p() {
    if(param.SM == 1) return 0.;

    double alphas_MSOFT = param.run.runningAlphasCalculation(param.MSOFT_Q);
    int nb_neut = (param.mass_neut[5] == 0.) ? 4 : 5;

    double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       (param.mu_Q + param.A_t / param.tan_beta) / param.mass_gluino *
                       (param.stop_mix[1][1] * param.stop_mix[1][1] * 
                        H2(param.mass_t2 * param.mass_t2 / param.mass_gluino / param.mass_gluino, param.mass_stl * param.mass_stl / param.mass_gluino / param.mass_gluino) +
                        param.stop_mix[1][2] * param.stop_mix[1][2] * 
                        H2(param.mass_t1 * param.mass_t1 / param.mass_gluino / param.mass_gluino, param.mass_stl * param.mass_stl / param.mass_gluino / param.mass_gluino));

    for(int ie = 1; ie <= nb_neut; ++ie) {
        epsilon0p += param.yub[3] * param.yub[3] / 16.0 / M_PI / M_PI * 
                     param.neut_mix[ie][4] * param.neut_mix[ie][3] * 
                     (param.mu_Q / param.tan_beta) / param.mass_neut[ie] *
                     (param.stop_mix[1][1] * param.stop_mix[1][1] * param.sbot_mix[1][1] * param.sbot_mix[1][1] * 
                      H2(param.mass_t1 * param.mass_t1 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b2 * param.mass_b2 / param.mass_neut[ie] / param.mass_neut[ie]) +
                      param.stop_mix[1][1] * param.stop_mix[1][1] * param.sbot_mix[1][2] * param.sbot_mix[1][2] * 
                      H2(param.mass_t1 * param.mass_t1 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b1 * param.mass_b1 / param.mass_neut[ie] / param.mass_neut[ie]) +
                      param.stop_mix[1][2] * param.stop_mix[1][2] * param.sbot_mix[1][1] * param.sbot_mix[1][1] * 
                      H2(param.mass_t2 * param.mass_t2 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b2 * param.mass_b2 / param.mass_neut[ie] / param.mass_neut[ie]) +
                      param.stop_mix[1][2] * param.stop_mix[1][2] * param.sbot_mix[1][2] * param.sbot_mix[1][2] * 
                      H2(param.mass_t2 * param.mass_t2 / param.mass_neut[ie] / param.mass_neut[ie], param.mass_b1 * param.mass_b1 / param.mass_neut[ie] / param.mass_neut[ie]));

    }

    return epsilon0p;
}


// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_1p
double EpsilonCalculator::epsilon_1p() const {
    if(param.SM == 1) return 0.;

    // Calcul du premier terme en utilisant yub[3], A_b, MqL3_Q, MbR_Q, mu_Q
    double term1 = 1.0 / 16.0 / M_PI / M_PI * 
                   (param.yub[3] * param.yub[3] * param.A_b / param.mu_Q * 
                    H2(std::pow(param.MqL3_Q / param.mu_Q, 2), std::pow(param.MbR_Q / param.mu_Q, 2)));

    // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
    double term2 = -param.g2 * param.g2 * param.M2_Q / param.mu_Q * 
                   H2(std::pow(param.MqL3_Q / param.mu_Q, 2), std::pow(param.M2_Q / param.mu_Q, 2)) / 16.0 / M_PI / M_PI;

    return term1 + term2;
}
