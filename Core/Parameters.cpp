#include "Parameters.h"
#include <iostream>

Parameters::Parameters() {
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
    masses[24] = 8.04229965E+01;
    masses[25] = 1.15104301E+02;
    coupling[1] = 3.57522130E-01;
    coupling[2] = 6.52355075E-01;

    extpar[25] = 10.; // tanb
    extpar[11] = -3800.; //At(MX)
    extpar[12] = -3800.; //Ab(MX)
    A_t = 0; 
    MqL3_Q = 0;
    MbR_Q = 0;
    mass_stl = 0;

    run = QCDParameters(masses[5], masses[5], masses[6], masses[6]);
}

Parameters *Parameters::GetInstance(int index)
{

    if(index < 0 || index > 1) {
            std::cerr << "Invalid Index. Must be 0 or 1." << std::endl;
            return nullptr;
        }

    if (!Parameters::instance[index]) {
        Parameters::instance[index] = new Parameters();
    }
    return Parameters::instance[index];
}

void Parameters::setScale(double Q) {

    // this->Q = Q;
    // this->sm.mass_b_Q = run.runningAlphasCalculation(Q);
    // this->sm.mass_t_Q = run.runningAlphasCalculation(Q);
}

Parameters* Parameters::instance[2] = {nullptr, nullptr};

