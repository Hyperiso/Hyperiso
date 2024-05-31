#include "Parameters.h"
#include "Logger.h"
#include "MemoryManager.h"
#include "Interface.h"
#include <iostream>
#include <complex>
#include <span>

std::string doubleToString(double value, int precision) {

	std::ostringstream out;
	out << std::fixed << std::setprecision(precision) << value;
	return out.str();
}

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
        std::vector<std::string> ids = matrixIds(SIZE);
        std::vector<double> values (SIZE * SIZE);
        lha->extractFromBlock(blockName, values, ids);
        for (int i=0; i!=values.size(); ++i) {
            matrix[ids[i][0] - 49][ids[i][2] - 49] = values[i];
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
        case 2:
            initTHDM();
            break;
        case 3: // Flavor
            initFlavor();
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

    sminputs[0] = 0.;
    sminputs[1] = *sm_inputs[0];
    Logger::getInstance()->debug("sminputs" + std::to_string(sminputs[1]));
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
    masses[6] = QCDRunner.get_mt_mt();       // t (running, m_t)
    masses[11] = 0.511e-3;      // e (pole)
    masses[13] = 0.105658;      // mu (pole)
    masses[15] = m_tau_pole;    // tau (pole)
    masses[23] = m_Z_pole;      // Z (pole)
    masses[24] = m_W;           // W  (running MW_MZ)
    masses[25] = 125.1;         // h0 

    Logger::getInstance()->info("mW : " + std::to_string(m_W));
    
    // Couplings
    double sW = std::sqrt(1 - std::pow(m_W / m_Z_pole, 2));

    std::cout << asin(sW) << std::endl;

    gauge[2] = std::pow(2, 1.25) * m_W * std::sqrt(G_F);     // g2
    gauge[1] = gauge[2] * sW / std::sqrt(1 - sW * sW);       // gp 
    gauge[3] = std::sqrt(4 * M_PI * alpha_s_MZ);             // gs
    gauge[4] = std::sqrt(4 * M_PI / inv_alpha_em);           // e_em     
    Logger::getInstance()->info("gp : " + std::to_string(gauge[1]));
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

    LhaReader* lha = MemoryManager::GetInstance()->getReader();

    std::string root = MemoryManager::findNearestHyperisoDirectory();
    std::string spectrumFile = root + "Test/spectrum.slha";
    Logger::getInstance()->info("Starting SUSY spectrum calculation...");
    
    CalculatorType calculatorType = CalculatorType::Softsusy;
    GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", lha->getLhaPath(), spectrumFile);

    lha->update(spectrumFile);
    Logger::getInstance()->info("LHA Blocks updated.");

    std::vector<std::string> mandatory {"STOPMIX", "SBOTMIX", "STAUMIX", "UMIX", "VMIX", "NMIX", "YU", "YD", "YE", "AU", "AD", "AE", "ALPHA", "HMIX", "GAUGE", "MSOFT", "MASS"};
    if (this->checkLHA(mandatory)) {
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

        this->alpha = lha->getValue<double>("ALPHA", "");
        this->susy_Q = static_cast<LhaElement<double>*>(lha->getBlock("HMIX")->get("1"))->getScale();
        
        std::vector<double> values (4);
        lha->extractFromBlock("HMIX", values);
        for (int i=0; i!=values.size(); ++i) {
            this->hmix[i + 1] = values[i];
        }
        
        values.resize(3);
        lha->extractFromBlock("GAUGE", values);
        for (int i=0; i!=values.size(); ++i) {
            this->gauge[i + 1] = values[i];
        }

        auto elts = lha->getBlock("MSOFT")->getEntries();
        for (size_t i = 0; i < elts->size(); ++i) {
            auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
            this->msoft[std::stoi(e->getId())] = e->getValue();
        }
        
        elts = lha->getBlock("MASS")->getEntries();
        for (size_t i = 0; i < elts->size(); ++i) {
            auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
            this->masses[std::stoi(e->getId())] = e->getValue();
        }
    } else {
        Logger::getInstance()->error("Cannot intialize SUSY parameters: LHA file is incomplete.");
    }

    
}

void Parameters::initFlavor() {
    // Hardcoded for now, should read from LHA file eventually
    this->masses[511] = 5.27958;
    this->masses[531] = 5.36677;

    this->lifetimes["511"] = 1.519e-12;
    this->lifetimes["531"] = 1.510e-12;

    this->fconst["511|1"] = 0.1905;
    this->fconst["531|1"] = 0.2277;
}

bool Parameters::checkLHA(std::vector<std::string> mandatory_blocks)
{
    LhaReader* lha = MemoryManager::GetInstance()->getReader();
    for (auto b: mandatory_blocks) {
        if (!lha->hasBlock(b))
            return false;
    }
    return true;
}

void Parameters::initTHDM() {
    LhaReader* lha = MemoryManager::GetInstance()->getReader();
    MemoryManager * memo = MemoryManager::GetInstance();

    std::string root = MemoryManager::findNearestHyperisoDirectory();
    std::string spectrumFile = root + "Test/thdm_spectrum.lha";
    Logger::getInstance()->info("Starting THDM spectrum calculation...");
    // TwoHDMCalculatorFactory::executeCommand("calculateSpectrum", memo->getInputLhaPath(), spectrumFile);

    CalculatorType calculatorType = CalculatorType::TwoHDM;
    GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", memo->getInputLhaPath(), spectrumFile);

    Logger::getInstance()->info("WAOUW : " + memo->getInputLhaPath());  
    
    lha->update(spectrumFile);
    
    Logger::getInstance()->info("LHA Blocks updated.");

    std::vector<std::string> mandatory {"MINPAR", "MASS", "ALPHA"};
    if (this->checkLHA(mandatory)) {
        // Read MASS block
        auto elts = lha->getBlock("MASS")->getEntries();
        for (size_t i = 0; i < elts->size(); ++i) {
            auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
            this->masses[std::stoi(e->getId())] = e->getValue();
        }
        
        // Read ALPHA block
        this->alpha = lha->getValue<double>("ALPHA", "");

        // Read tan(beta) from block MINPAR
        double tan_beta = lha->getValue<double>("MINPAR", "3");
        this->hmix[2] = tan_beta;

        // Read yukawa type and populate the yukawa blocks appropriately
        int type = static_cast<int>(lha->getValue<double>("MINPAR", "24"));
        double cot_beta = 1 / tan_beta;
        this->yu[2][2] = cot_beta;
        switch (type) {
            case 1:
                this->yd[2][2] = cot_beta;
                this->ye[2][2] = cot_beta;
                break;
            case 2:
                this->yd[2][2] = -tan_beta;
                this->ye[2][2] = -tan_beta;
                break;
            case 3:
                this->yd[2][2] = -tan_beta;
                this->ye[2][2] = cot_beta;
                break;
            case 4:
                this->yd[2][2] = cot_beta;
                this->ye[2][2] = -tan_beta;
                break;
            default:
                Logger::getInstance()->error("Cannot initialize THDM parameters: Unknown Yukawa type " + std::to_string(type));
        }
        Logger::getInstance()->info("THDM parameters initialized.");
    } else {
        Logger::getInstance()->error("Cannot intialize THDM parameters: LHA file is incomplete.");
    }
}

Parameters *Parameters::GetInstance(int modelId)
{

    if(modelId < 0 || modelId > 3) {
        std::cerr << "Invalid Index. Must be 0, 1 or 2." << std::endl;
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

double Parameters::running_mass(double quarkmass, double Q_init, double Q_end,  std::string option_massb, std::string option_masst) {
    if (quarkmass >= masses[4]) {
        return this->QCDRunner.running_mass(quarkmass, Q_init, Q_end, option_massb, option_masst);
    } else {
        Logger::getInstance()->error("In Parameters::running_mass: Quark of mass " + std::to_string(quarkmass) + " lower than charm mass, not possible.");
        return 0;
    }
}

double return_if_defined(std::map<std::string, double>& map, const std::string& id, const std::string& error_label) {
    if (map.contains(id)) {
        return map[id];
    } else {
        Logger::getInstance()->warn(error_label + " with key [" + id + "] is undefined.");
        return NAN;
    }
}

double Parameters::getFlavorParam(FlavorParamType type, const std::string& id) {
    switch (type) {
        case FlavorParamType::LIFETIME:
            return return_if_defined(this->lifetimes, id, "Lifetime");
        case FlavorParamType::DECAY_CONSTANT:
            return return_if_defined(this->fconst, id, "Decay constant");
        default:
            Logger::getInstance()->error("Unknown parameter type.");
            return NAN;
    }
}

Parameters* Parameters::instance[4] = {nullptr, nullptr, nullptr, nullptr};
