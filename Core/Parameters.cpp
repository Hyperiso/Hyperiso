#include "Parameters.h"

Parameters::Parameters() {
    // Initialisation des paramètres utilisés dans les calculs epsilon_x
    sm.SM = 1;  // Exemple de valeur, à ajuster selon votre cas
    sm.gp = 0;
    sm.g2 = 0;
    sm.MSOFT_Q = 0;
    sm.mass_top_pole = 0;
    sm.mass_b_pole = 0;
    sm.mass_b_Q = 0;
    sm.mass_t_Q = 0;
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

    masses[6] = 173.0;  // Top quark
    masses[5] = 4.18;
    coupling[1] = 1;
    coupling[2] = 1;
    A_t = 0;
    MqL3_Q = 0;
    MbR_Q = 0;
    mass_stl = 0;
}

Parameters *Parameters::GetInstance()
{
    if (!Parameters::instance) {
        Parameters::instance = new Parameters();
    }
    return Parameters::instance;
}

void Parameters::setScale(double Q) {

    this->Q = Q;
    this->sm.mass_b_Q = run.runningAlphasCalculation(Q);
    this->sm.mass_t_Q = run.runningAlphasCalculation(Q);
}
