#include "Parameters.h"
#include "Logger.h"
#include "MemoryManager.h"
#include <iostream>
#include <complex>

typedef std::complex<double> complex_t; 

Parameters::Parameters(int modelId) {
    switch (modelId) {
        case 0: // SM
            initSM();
            break;
        case 1: // SUSY
            initSUSY();
            break;
        default:
            Logger::getInstance()->error("Trying to instantiate parameters for unknown model ID " + std::to_string(modelId));
    }
}

void Parameters::initSM() {

    /* 
        Reading LHA blocks if given 
    */

    LhaReader* lha = MemoryManager::GetInstance()->getReader();

    // SMINPUTS
    double inv_alpha_em{1.27934e2}, G_F{1.16637e-5}, alpha_s_MZ{1.184e-1}, m_Z_pole{91.1876}, m_b_mb{4.18}, m_t_pole{172.7}, m_tau_pole{1.777};
    std::vector<double*> sm_inputs = {&inv_alpha_em, &G_F, &alpha_s_MZ, &m_Z_pole, &m_b_mb, &m_t_pole, &m_tau_pole};
    lha->extractFromBlock("SMINPUTS", sm_inputs);

    // VCKMIN 
    double lambda{0.22500}, A{0.826}, rho{0.159}, eta{0.348};
    std::vector<double*> ckm_inputs = {&lambda, &A, &rho, &eta};
    lha->extractFromBlock("VCKMIN", ckm_inputs);

    // TODO: PMNS matrix

    /* 
        Initializing SM parameters
    */

    double m_W = std::sqrt(std::pow(m_Z_pole, 2) / 2 + std::sqrt(std::pow(m_Z_pole, 4) / 4 - M_PI * std::pow(m_Z_pole, 2) / inv_alpha_em / G_F / std::sqrt(2)));
    QCDRunner = QCDParameters(alpha_s_MZ, m_Z_pole, m_b_mb, m_b_mb, m_t_pole, m_t_pole);

    // Masses (from PDG 2023)
    masses[1] = 4.7e-3;         // d (2 GeV)
    masses[2] = 2.2e-3;         // u (2 GeV)
    masses[3] = 93e-3;          // s (2 GeV)
    masses[4] = 1.27;           // c (running, 1.27 GeV)
    masses[5] = 4.18;           // b (running, 4.18 GeV)
    masses[6] = m_t_pole;       // t (pole)
    masses[11] = 0.511e-3;      // e (pole)
    masses[13] = 0.105658;      // mu (pole)
    masses[15] = m_tau_pole;    // tau (pole)
    masses[23] = m_Z_pole;      // Z (running MZ_MZ)
    masses[24] = m_W;           // W  (running MW_MZ)
    masses[25] = 125.1;         // h0 

    // Couplings
    double sW = std::sqrt(1 - std::pow(m_W / m_Z_pole, 2));
    coupling[2] = std::pow(2, 1.25) * m_W * std::sqrt(G_F);   // g2
    coupling[1] = coupling[2] * sW / std::sqrt(1 - sW * sW);   // gp 
    coupling[3] = std::sqrt(4 * M_PI * alpha_s_MZ); // gs
    coupling[4] = std::sqrt(4 * M_PI / inv_alpha_em); // e_em     

    // CKM Matrix
    ckm[0][0] = 1 - lambda * lambda / 2;
    ckm[0][1] = lambda;
    ckm[0][2] = A * lambda * lambda * lambda * complex_t{rho, -eta};
    ckm[1][0] = -lambda;
    ckm[1][1] = 1 - lambda * lambda / 2;
    ckm[1][2] = A * lambda * lambda;
    ckm[2][0] = A * lambda * lambda * lambda * complex_t{1 - rho, -eta};
    ckm[2][1] = -A * lambda * lambda;
    ckm[2][2] = 1;

}

void Parameters::initSUSY() {
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

    extpar[25] = 10.; // tanb
    extpar[11] = -3800.; //At(MX)
    extpar[12] = -3800.; //Ab(MX)
    A_t = 0; 
    MqL3_Q = 0;
    MbR_Q = 0;
    mass_stl = 0;
}

Parameters *Parameters::GetInstance(int modelId)
{

    if(modelId < 0 || modelId > 1) {
        std::cerr << "Invalid Index. Must be 0 or 1." << std::endl;
        return nullptr;
    }

    if (!Parameters::instance[modelId]) {
        Parameters::instance[modelId] = new Parameters(modelId);
    }
    return Parameters::instance[modelId];
}

void Parameters::setScale(double Q) {

    // this->Q = Q;
    // this->sm.mass_b_Q = run.runningAlphasCalculation(Q);
    // this->sm.mass_t_Q = run.runningAlphasCalculation(Q);
}

Parameters* Parameters::instance[2] = {nullptr, nullptr};

