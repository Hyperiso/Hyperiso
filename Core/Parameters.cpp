#include "Parameters.h"

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
        for (size_t i=0; i!=values.size(); ++i) {
            matrix[ids[i][0] - 49][ids[i][2] - 49] = values[i];
        }
    } 
}

std::map<int, Parameters*> Parameters::instances;
std::map<int, Parameters*> ParametersFactory::instances;

Parameters* Parameters::GetInstance(int modelId) {
    return ParametersFactory::GetParameters(modelId);
}

Parameters::Parameters(ModelStrategy* modelStrategy)
    : strategy(modelStrategy) {
    strategy->initializeParameters(*this);
}

double Parameters::operator()(const std::string& block, int pdgCode) {
    return blockAccessor.getValue(block, pdgCode);
}

double Parameters::alpha_s(double Q) {
    return this->QCDRunner.runningAlphasCalculation(Q);
}

double Parameters::running_mass(double quarkmass, double Q_init, double Q_end,  std::string option_massb, std::string option_masst) {
    return this->QCDRunner.running_mass(quarkmass, Q_init, Q_end, option_massb, option_masst);
}

double Parameters::get_QCD_masse(std::string masstype) {
    if (masstype == "mt_mt"){
        return this->QCDRunner.get_mt_mt();
    }
    if (masstype == "mb_pole") {
        return this->QCDRunner.get_mb_pole();
    }
    if (masstype == "mt_pole") {
        return this->QCDRunner.get_mt_pole();
    }
    if (masstype == "mb_mb") {
        return this->QCDRunner.get_mb_mb();
    }
    if (masstype == "mb_1S") {
        return this->QCDRunner.mb_1S();
    }
}

// SMModelStrategy implementation
void SMModelStrategy::initializeParameters(Parameters& params) {

    LhaReader* lha = MemoryManager::GetInstance()->getReader();

    // SMINPUTS
    double inv_alpha_em{1.37934e2}, G_F{1.16637e-5}, alpha_s_MZ{1.184e-1}, m_Z_pole{91.1876}, m_b_mb{4.18}, m_t_pole{172.9}, m_tau_pole{1.777};
    std::vector<double*> sm_inputs = {&inv_alpha_em, &G_F, &alpha_s_MZ, &m_Z_pole, &m_b_mb, &m_t_pole, &m_tau_pole};
    lha->extractFromBlock("SMINPUTS", sm_inputs);
    params.addBlock("SMINPUTS", std::make_unique<SMInputBlock>());

    params.setBlockValue("SMINPUTS", 0, 0.);
    params.setBlockValue("SMINPUTS", 1, *sm_inputs[0]);
    params.setBlockValue("SMINPUTS", 2, *sm_inputs[1]);
    params.setBlockValue("SMINPUTS", 5, *sm_inputs[4]);

    // VCKMIN 
    params.addBlock("RECKM", std::make_unique<CKMBlock>());
    params.addBlock("IMCKM", std::make_unique<CKMBlock>());
    double lambda{0.22500}, A{0.826}, rho{0.159}, eta{0.348};
    std::vector<double*> ckm_inputs = {&lambda, &A, &rho, &eta};
    lha->extractFromBlock("VCKMIN", ckm_inputs);
    params.setBlockValue("RECKM", 0, std::real(1 - lambda * lambda / 2));
    params.setBlockValue("IMCKM", 0, std::imag(1 - lambda * lambda / 2));

    params.setBlockValue("RECKM", 1, std::real(lambda));
    params.setBlockValue("IMCKM", 1, std::imag(lambda));

    params.setBlockValue("RECKM", 2, std::real(A * lambda * lambda * lambda * complex_t{rho, -eta}));
    params.setBlockValue("IMCKM", 2, std::imag(A * lambda * lambda * lambda * complex_t{rho, -eta}));

    params.setBlockValue("RECKM", 10, std::real(-lambda));
    params.setBlockValue("IMCKM", 10, std::imag(-lambda));

    params.setBlockValue("RECKM", 11, std::real(1 - lambda * lambda / 2));
    params.setBlockValue("IMCKM", 11, std::imag(1 - lambda * lambda / 2));

    params.setBlockValue("RECKM", 12, std::real(A * lambda * lambda));
    params.setBlockValue("IMCKM", 12, std::imag(A * lambda * lambda));

    params.setBlockValue("RECKM", 20, std::real(A * lambda * lambda * lambda * complex_t{1 - rho, -eta}));
    params.setBlockValue("IMCKM", 20, std::imag(A * lambda * lambda * lambda * complex_t{1 - rho, -eta}));

    params.setBlockValue("RECKM", 21, std::real(-A * lambda * lambda));
    params.setBlockValue("IMCKM", 21, std::imag(-A * lambda * lambda));

    params.setBlockValue("RECKM", 22, std::real(1));
    params.setBlockValue("IMCKM", 22, std::imag(1));


    double m_W = std::sqrt(std::pow(m_Z_pole, 2) / 2 + std::sqrt(std::pow(m_Z_pole, 4) / 4 - M_PI * std::pow(m_Z_pole, 2) / inv_alpha_em / G_F / std::sqrt(2)));
    //Masses
    params.addBlock("MASS", std::make_unique<MassBlock>());
    // Masses (from PDG 2023)
    params.setBlockValue("MASS", 1, 4.7e-3);
    params.setBlockValue("MASS", 2, 2.2e-3);
    params.setBlockValue("MASS", 3, 93e-3);
    params.setBlockValue("MASS", 4, 1.27);

    params.setQCDParameters(QCDParameters(alpha_s_MZ, m_Z_pole, m_t_pole, m_b_mb, params("MASS", 2), params("MASS", 1), params("MASS", 3), params("MASS", 4)));

    params.setBlockValue("MASS", 5, m_b_mb);
    params.setBlockValue("MASS", 6, params.get_QCD_masse("mt_mt"));
    params.setBlockValue("MASS", 11, 0.511e-3);
    params.setBlockValue("MASS", 13, 0.105658);
    params.setBlockValue("MASS", 15, m_tau_pole);
    params.setBlockValue("MASS", 23, m_Z_pole);
    params.setBlockValue("MASS", 24, m_W);
    params.setBlockValue("MASS", 25, 125.1);


    // Couplings
    double sW = std::sqrt(1 - std::pow(m_W / m_Z_pole, 2));
    params.addBlock("GAUGE", std::make_unique<GaugeBlock>());

    params.setBlockValue("GAUGE", 2, std::pow(2, 1.25) * m_W * std::sqrt(G_F));  // g2
    params.setBlockValue("GAUGE", 1, params("GAUGE", 2) * sW / std::sqrt(1 - sW * sW)); // gp 
    params.setBlockValue("GAUGE", 3, std::sqrt(4 * M_PI * alpha_s_MZ)); // gs
    params.setBlockValue("GAUGE", 4, std::sqrt(4 * M_PI / inv_alpha_em)); // e_em 

}

// SUSYModelStrategy implementation
void SUSYModelStrategy::initializeParameters(Parameters& params) {

    LhaReader* lha = MemoryManager::GetInstance()->getReader();

    std::string root = project_root.data();
    std::string spectrumFile = root +"/" + "Test/spectrum.slha";
    LOG_INFO("Starting SUSY spectrum calculation...");
    
    CalculatorType calculatorType = CalculatorType::Softsusy;
    GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", lha->getLhaPath(), spectrumFile);
    
    LOG_INFO("SUSY spectrum calculation ran sucessfully");

    lha->update(spectrumFile);
    LOG_INFO("LHA Blocks updated.");
    

    std::vector<std::string> mandatory {"STOPMIX", "SBOTMIX", "STAUMIX", "UMIX", "VMIX", "NMIX", "YU", "YD", "YE", "AU", "AD", "AE", "ALPHA", "HMIX", "GAUGE", "MSOFT", "MASS"};

    std::map<std::string, std::shared_ptr<ArrayBlock<2,2>>> type22 {
        {"STOPMIX", std::make_shared<StopMixBlock>()},
        {"SBOTMIX", std::make_shared<SbotMixBlock>()},
        {"STAUMIX", std::make_shared<StauMixBlock>()},
        {"UMIX", std::make_shared<UMIXBlock>()},
        {"VMIX", std::make_shared<VMIXBlock>()}
    };

    std::map<std::string, std::shared_ptr<ArrayBlock<4,4>>> type44 {
        {"NMIX", std::make_shared<NMIXBlock>()},
        {"H0MIX", std::make_shared<H0mixBlock>()},
        {"A0MIX", std::make_shared<A0mixBlock>()}
    };

    std::map<std::string, std::shared_ptr<ArrayBlock<3,3>>> type33 {
        {"YU", std::make_shared<YUBlock>()},
        {"YD", std::make_shared<YDBlock>()},
        {"YE", std::make_shared<YEBlock>()},
        {"AU", std::make_shared<AUBlock>()},
        {"AD", std::make_shared<ADBlock>()},
        {"AE", std::make_shared<AEBlock>()}
    };



    for (auto& elem : type22) {
        std::array<std::array<double, 2>,2> temp;
        readMatrix(temp, elem.first, lha);
        *(elem.second) = temp;
        std::unique_ptr<ArrayBlock<2, 2>> uniquePtr = std::make_unique<ArrayBlock<2,2>>(std::move(*elem.second));
        params.addBlock(elem.first, std::move(uniquePtr));
    }

    for (auto& elem : type44) {
        std::array<std::array<double, 4>,4> temp;
        readMatrix(temp, elem.first, lha);
        *(elem.second) = temp;
        std::unique_ptr<ArrayBlock<4, 4>> uniquePtr = std::make_unique<ArrayBlock<4,4>>(std::move(*elem.second));
        params.addBlock(elem.first, std::move(uniquePtr));
    }

    for (auto& elem : type33) {
        std::array<std::array<double, 3>,3> temp;
        readMatrix(temp, elem.first, lha);
        *(elem.second) = temp;
        std::unique_ptr<ArrayBlock<3, 3>> uniquePtr = std::make_unique<ArrayBlock<3,3>>(std::move(*elem.second));
        params.addBlock(elem.first, std::move(uniquePtr));
    }

    auto alphablock = std::make_unique<AlphaBlock>();
    alphablock->setValue(0, lha->getValue<double>("ALPHA", ""));
    params.addBlock("ALPHA", std::move(alphablock));

    auto hmixblock = std::make_unique<HMIXBlock>();
    hmixblock->setValue(0, static_cast<LhaElement<double>*>(lha->getBlock("HMIX")->get("1"))->getScale());
    std::vector<double> values (4);
    lha->extractFromBlock("HMIX", values);
    for (size_t i=0; i!=values.size(); ++i) {
        hmixblock->setValue(i+1, values[i]);
    }
    params.addBlock("HMIX", std::move(hmixblock));

    values.resize(3);
    auto gaugeblock = std::make_unique<GaugeBlock>();
    lha->extractFromBlock("GAUGE", values);
    for (size_t i=0; i!=values.size(); ++i) {
        gaugeblock->setValue(i+1, values[i]);
    }
    params.addBlock("GAUGE", std::move(gaugeblock));

    auto msoftblock = std::make_unique<MSOFTBlock>();
    auto elts = lha->getBlock("MSOFT")->getEntries();
    for (size_t i = 0; i < elts->size(); ++i) {
        auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
        msoftblock->setValue(std::stoi(e->getId()), e->getValue());
    }
    params.addBlock("MSOFT", std::move(msoftblock));

    auto massblock = std::make_unique<MassBlock>();
    elts = lha->getBlock("MASS")->getEntries();
    massblock->setValue(1000039, 0.);
    massblock->setValue(45, 0.);
    massblock->setValue(46, 0.);
    for (size_t i = 0; i < elts->size(); ++i) {
        auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
        massblock->setValue(std::stoi(e->getId()), e->getValue());
    }
    params.addBlock("MASS", std::move(massblock));

}

// THDMModelStrategy implementation
void THDMModelStrategy::initializeParameters(Parameters& params) {
    LhaReader* lha = MemoryManager::GetInstance()->getReader();
    MemoryManager * memo = MemoryManager::GetInstance();

    std::string root = project_root.data();
    std::string spectrumFile = root + "/" + "Test/thdm_spectrum.lha";
    LOG_INFO("Starting THDM spectrum calculation...");
    // TwoHDMCalculatorFactory::executeCommand("calculateSpectrum", memo->getInputLhaPath(), spectrumFile);

    CalculatorType calculatorType = CalculatorType::TwoHDM;
    GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", memo->getInputLhaPath(), spectrumFile);
    
    LOG_INFO("THDM spectrum calculation ran sucessfully");
 
    
    lha->update(spectrumFile);
    
    LOG_INFO("LHA Blocks updated.");

    std::vector<std::string> mandatory {"MINPAR", "MASS", "ALPHA"};

    auto massblock = std::make_unique<MassBlock>();
    auto elts = lha->getBlock("MASS")->getEntries();
    for (size_t i = 0; i < elts->size(); ++i) {
        auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
        massblock->setValue(std::stoi(e->getId()), e->getValue());
        // this->masses[std::stoi(e->getId())] = e->getValue();
    }
    params.addBlock("MASS", std::move(massblock));

    auto alphablock = std::make_unique<AlphaBlock>();
    alphablock->setValue(0, lha->getValue<double>("ALPHA", ""));
    params.addBlock("ALPHA", std::move(alphablock));

    auto hmixblock = std::make_unique<HMIXBlock>();
    double tan_beta = lha->getValue<double>("MINPAR", "3");
    hmixblock->setValue(2, tan_beta);
    params.addBlock("HMIX", std::move(hmixblock));

    int type = static_cast<int>(lha->getValue<double>("MINPAR", "24"));
    double cot_beta = 1 / tan_beta;

    auto yublock = std::make_unique<YUBlock>();
    auto ydblock = std::make_unique<YUBlock>();
    auto yeblock = std::make_unique<YUBlock>();

    yublock->setValue(22, cot_beta);

    
    switch (type) {
        case 1:
            ydblock->setValue(22, cot_beta);
            yeblock->setValue(22, cot_beta);
            yeblock->setValue(11, cot_beta);
            yeblock->setValue(00, cot_beta);
            break;
        case 2:
            ydblock->setValue(22, -tan_beta);
            yeblock->setValue(22, -tan_beta);
            yeblock->setValue(00, -tan_beta);
            yeblock->setValue(00, -tan_beta);
            break;
        case 3:
            ydblock->setValue(22, -tan_beta); 
            yeblock->setValue(22, cot_beta);
            yeblock->setValue(11, cot_beta);
            yeblock->setValue(00, cot_beta);
            break;
        case 4:
            ydblock->setValue(22, cot_beta);
            yeblock->setValue(22, -tan_beta);
            yeblock->setValue(11, -tan_beta);
            yeblock->setValue(00, -tan_beta);
            break;
        default:
            LOG_ERROR("ValueError", "Cannot initialize THDM parameters: Unknown Yukawa type " + std::to_string(type));
    }

    params.addBlock("YU", std::move(yublock));
    params.addBlock("YD", std::move(ydblock));
    params.addBlock("YL", std::move(yeblock));

}

// FlAVORModelStrategy implementation
void FlAVORModelStrategy::initializeParameters(Parameters& params) {

}

double return_if_defined(std::map<std::string, double>& map, const std::string& id, const std::string& error_label) {
    if (map.contains(id)) {
        return map[id];
    } else {
        LOG_WARN(error_label + " with key [" + id + "] is undefined.");
        return NAN;
    }
}

double Parameters::getFlavorParam(FlavorParamType type, const std::string& id) {
    return this->flavorblockAccessor.getValue(type, id);
    // switch (type) {

    //     case FlavorParamType::LIFETIME:
    //         return this->flavorblockAccessor.getValue(type, id);
    //         return return_if_defined(this->lifetimes, id, "Lifetime");
    //     case FlavorParamType::DECAY_CONSTANT:
    //         return return_if_defined(this->fconst, id, "Decay constant");
    //     default:
    //         LOG_ERROR("ValueError", "Unknown parameter type.");
    //         return NAN;
    // }
}

Parameters* ParametersFactory::GetParameters(int modelId) {
    if (instances.find(modelId) == instances.end()) {
        ModelStrategy* strategy = createStrategy(modelId);
        instances[modelId] = new Parameters(strategy);
    }
    return instances[modelId];
}

ModelStrategy* ParametersFactory::createStrategy(int modelId) {
    switch (modelId) {
        case 0:
            return new SMModelStrategy();
        case 1:
            return new SUSYModelStrategy();
        case 2:
            return new THDMModelStrategy();
        case 3:
            return new FlAVORModelStrategy();
        // Add other cases for different models (THDM, Flavor, etc.)
        default:
            throw std::invalid_argument("Unknown model ID");
    }
}
