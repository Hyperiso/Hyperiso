#include "Parameters.h"
#include "Logger.h"
#include "MemoryManager.h"
#include <iostream>
#include <complex>
#include <span>

typedef std::complex<double> complex_t; 

std::vector<std::string> matrixIds(int size) {
    std::vector<std::string> ids;
    for (int i = 1; i <= size; ++i) {
        for (int j = 1; j <= size; ++j) {
            std::stringstream ss;
            ss << i << "|" << j;
            ids.emplace_back(ss.str());
        }
    }
    return ids;
}

template<std::size_t SIZE>
void readMatrix(std::array<std::array<double, SIZE>, SIZE>& matrix, std::string blockName, LhaReader* lha) {
    if (lha->hasBlock(blockName)) {
        auto ids = matrixIds(SIZE);
        std::vector<double*> values;
        values.resize(SIZE * SIZE);
        lha->extractFromBlock(blockName, values, ids);
        for (int i=0; i!=values.size(); ++i) {
            matrix[30 + ids[i][0]][30 + ids[i][2]] = values[i] ? *values[i] : 0;
        }
    } 
}

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
    double inv_alpha_em{1.37934e2}, G_F{1.16637e-5}, alpha_s_MZ{1.184e-1}, m_Z_pole{91.1876}, m_b_mb{4.18}, m_t_pole{172.9}, m_tau_pole{1.777};
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

    // Masses (from PDG 2023)
    masses[1] = 4.7e-3;         // d (2 GeV)
    masses[2] = 2.2e-3;         // u (2 GeV)
    masses[3] = 93e-3;          // s (2 GeV)
    masses[4] = 1.27;           // c (running, 1.27 GeV)

    QCDRunner = QCDParameters(alpha_s_MZ, m_Z_pole, m_t_pole, m_b_mb, masses[2],masses[1],masses[3],masses[4]);

    masses[5] = 4.18;           // b (running, 4.18 GeV)
    masses[6] = QCDRunner.get_mt_mt();       // t (pole)
    masses[11] = 0.511e-3;      // e (pole)
    masses[13] = 0.105658;      // mu (pole)
    masses[15] = m_tau_pole;    // tau (pole)
    masses[23] = m_Z_pole;      // Z (pole)
    masses[24] = m_W;           // W  (running MW_MZ)
    masses[25] = 125.1;         // h0 

    
    // Couplings
    double sW = std::sqrt(1 - std::pow(m_W / m_Z_pole, 2));

    gauge[2] = std::pow(2, 1.25) * m_W * std::sqrt(G_F);     // g2
    gauge[1] = gauge[2] * sW / std::sqrt(1 - sW * sW);    // gp 
    gauge[3] = std::sqrt(4 * M_PI * alpha_s_MZ);             // gs
    gauge[4] = std::sqrt(4 * M_PI / inv_alpha_em);           // e_em     

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

    /*
        This is where we should launch any spectrum calculation code to update the SLHA file with the corresponding blocks

        SpectrumCalculator.calculate(modsel, minpar, extpar);
    */

   LhaReader* lha = MemoryManager::GetInstance()->getReader();

    readMatrix(this->stopmix, "STOPMIX", lha);
    readMatrix(this->sbotmix, "SBOTMIX", lha);
    readMatrix(this->staumix, "STAUMIX", lha);
    readMatrix(this->umix, "UMIX", lha);
    readMatrix(this->vmix, "VMIX", lha);
    readMatrix(this->nmix, "NMIX", lha);
    readMatrix(this->yu, "YU", lha);
    readMatrix(this->yd, "YD", lha);
    readMatrix(this->ye, "YE", lha);
    readMatrix(this->au, "AU", lha);
    readMatrix(this->ad, "AD", lha);
    readMatrix(this->ae, "AE", lha);

    if (lha->hasBlock("ALPHA")) {
        std::vector<double*> values = {&(this->alpha)};
        lha->extractFromBlock("ALPHA", values);
    } 

    if (lha->hasBlock("HMIX")) {
        std::vector<double*> values;
        values.resize(4);
        lha->extractFromBlock("HMIX", values);
        for (int i=0; i!=values.size(); ++i) {
            this->hmix[i + 1] = *values[i];
        }

        this->susy_Q = static_cast<LhaElement<double>*>(lha->getBlock("HMIX")->get("1"))->getScale();
    } 

    if (lha->hasBlock("GAUGE")) {
        std::vector<double*> values;
        values.resize(3);
        lha->extractFromBlock("GAUGE", values);
        for (int i=0; i!=values.size(); ++i) {
            this->gauge[i + 1] = *values[i];
        }
    } 

    if (lha->hasBlock("MSOFT")) {
        auto elts = lha->getBlock("MSOFT")->getEntries();
        for (size_t i = 0; i < elts->size(); ++i) {
            auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
            this->msoft[std::stoi(e->getId())] = e->getValue();
        }
    } 
    
    if (lha->hasBlock("MASS")) {
        auto elts = lha->getBlock("MASS")->getEntries();
        for (size_t i = 0; i < elts->size(); ++i) {
            auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
            this->masses[std::stoi(e->getId())] = e->getValue();
        }
    } 
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

double Parameters::alpha_s(double Q) {
    return this->QCDRunner.runningAlphasCalculation(Q);
}

double Parameters::running_mass(double quark_m, double q_init, double q_fin, std::string opt_mb, std::string opt_mt) {
    return this->QCDRunner.running_mass(quark_m, q_init, q_fin, opt_mb, opt_mt);
}

Parameters* Parameters::instance[2] = {nullptr, nullptr};

