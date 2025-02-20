#include "Parameters.h"
#include <ranges>
#include <algorithm>

std::string doubleToString(double value, int precision) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(precision) << value;
	return out.str();
}

typedef std::complex<double> complex_t; 

std::vector<LhaID> matrixIds(int size) {
    std::vector<LhaID> ids;
    for (int i = 1; i <= size; ++i) {
        for (int j = 1; j <= size; ++j) {
            ids.emplace_back(LhaID({i, j}));
        }
    }
    return ids;
}

template<std::size_t SIZE>
void readMatrix(std::array<std::array<double, SIZE>, SIZE>& matrix, std::string blockName, LhaReader* lha) {
    if (lha->hasBlock(blockName)) {
        auto ids = matrixIds(SIZE);
        std::vector<double> values (SIZE * SIZE);
        lha->extractFromBlock(blockName, values, ids);
        for (size_t i=0; i!=values.size(); ++i) {
            matrix[ids[i].parts[0] - 49][ids[i].parts[1] - 49] = values[i];
        }
    } 
}

void populate_from_json(std::shared_ptr<MapBlock> block, std::vector<Value> json_vals) {
    for (const auto& val : json_vals) {
        if (val.name.starts_with(block->blockname)) {
            std::string del = "|";
            auto split = val.name.find(del);
            int pdg = std::stoi(val.name.substr(split + 1, val.name.size() - split));
            block->setValue(pdg, val.central_value);
        }
    }
}

void overwrite_from_lha(std::shared_ptr<MapBlock> block, LhaReader* reader) {
    if (!reader->hasBlock(block->blockname)) 
        return;

    for (const std::vector<int> sub_ids : LhaParamsHelper::get_minimal_content(block->blockname)) {
        LhaID id {sub_ids};
        if (reader->hasElement(block->blockname, id)) {
            block->setValue(id, reader->getValue<double>(block->blockname, id));
        }
    }
}

std::map<ParameterType, std::shared_ptr<Parameters>> Parameters::instances;
std::map<ParameterType, std::shared_ptr<Parameters>> ParametersFactory::instances;

std::shared_ptr<Parameters> Parameters::GetInstance(ParameterType id) {
    auto allowed = MemoryManager::GetInstance()->getParameterTypes();
    if (std::find(allowed.begin(), allowed.end(), id) == allowed.end())
        LOG_ERROR("OutOfRange", "Parameter type undefined");
    return ParametersFactory::GetParameters(id);
}

void Parameters::CleanupInstance(ParameterType id) {
    ParametersFactory::removeParameters(id);
}

ParameterType Parameters::GetType(const std::string &block, int pdgCode) {
    auto allowed_types = MemoryManager::GetInstance()->getParameterTypes();
    for (ParameterType tp : allowed_types) {
        auto p = Parameters::GetInstance(tp);
        for (auto &b : p->get_blocks_list()) {
            if (b == block) {
                if (p->get_block_infos(b).contains(pdgCode))
                    return tp;
            }
        }
    }
    LOG_ERROR("Invalid Parameter", "Parameter", block, ",", pdgCode, "is undefined.");
}

double Parameters::Get(ParameterType type, const std::string& block, int code) {
    LOG_VERBOSE("Attempting to retrieve parameter from instance", (int)type, "with id (", block, ",", code, ")");
    return (*GetInstance(type))(block, code);
}

double Parameters::Get(ParamId id) {
    return Get(id.type, id.block, id.code);
}

Parameters::Parameters(std::shared_ptr<ModelStrategy> modelStrategy)
    : strategy(modelStrategy) { 
    LOG_VERBOSE("Param creation at", this);
    strategy->initializeParameters(*this);
}

double Parameters::operator()(const std::string& block, int pdgCode) {
    return blockAccessor.getValue(block, pdgCode);
}

bool Parameters::exist(const std::string& block, int pdgCode) {
    return blockAccessor.exist(block, pdgCode);
}

void Parameters::addBlock(const std::string& name, std::shared_ptr<Block> block) {
    blockAccessor.addBlock(name, block);
}

void Parameters::setBlockValue(const std::string& name, int pdgCode, double value, bool force) {
    blockAccessor.setValue(name, pdgCode, value, force);
}

std::map<int, double> Parameters::get_block_infos(std::string blockName) {
    return blockAccessor.getAllValues(blockName);
}


std::vector<std::string> Parameters::get_blocks_list() {
    return blockAccessor.get_blocks();
}

complex_t Parameters::get_c_CKM_entry(int idx) {
    auto p = Parameters::GetInstance();
    return complex_t((*p)("RECKM", idx), (*p)("IMCKM", idx));
}

void SMModelStrategy::initializeParameters(Parameters& params) {

    LhaReader* lha = MemoryManager::GetInstance()->getReader();

    // SMINPUTS
    double inv_alpha_em{1.37934e2}, G_F{1.16637e-5}, alpha_s_MZ{1.184e-1}, m_Z_pole{91.1876}, m_b_mb{4.18}, m_t_pole{172.9}, m_tau_pole{1.777};
    std::vector<double*> sm_inputs = {&inv_alpha_em, &G_F, &alpha_s_MZ, &m_Z_pole, &m_b_mb, &m_t_pole, &m_tau_pole};
    lha->extractFromBlock("SMINPUTS", sm_inputs);
    params.addBlock("SMINPUTS", std::make_shared<SMInputBlock>());
    for (size_t i = 0; i < 7; i++) {
        params.setBlockValue("SMINPUTS", i + 1, *sm_inputs[i]);
    }
    params.setBlockValue("SMINPUTS", 10, 0.313);

    // VCKMIN 
    params.addBlock("RECKM", std::make_shared<RECKMBlock>());
    params.addBlock("IMCKM", std::make_shared<IMCKMBlock>());
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
    
    //Masses (from PDG 2023)
    params.addBlock("MASS", std::make_shared<MassBlock>());
    params.setBlockValue("MASS", 1, 4.7e-3);
    params.setBlockValue("MASS", 2, 2.2e-3);
    params.setBlockValue("MASS", 3, 93e-3);
    params.setBlockValue("MASS", 4, 1.27);

    QCDHelper::Init(alpha_s_MZ, m_Z_pole, m_t_pole, m_b_mb, 1.27, 93e-3, 2.2e-3, 4.7e-3);
    params.setBlockValue("MASS", 5, m_b_mb);
    params.setBlockValue("MASS", 6, m_t_pole);

    params.setBlockValue("MASS", 11, 0.511e-3);
    params.setBlockValue("MASS", 13, 0.105658);
    params.setBlockValue("MASS", 15, m_tau_pole);
    params.setBlockValue("MASS", 23, m_Z_pole);
    params.setBlockValue("MASS", 24, m_W);
    params.setBlockValue("MASS", 25, 125.1);

    // Couplings
    double sW = std::sqrt(1 - std::pow(m_W / m_Z_pole, 2));
    params.addBlock("GAUGE", std::make_shared<GaugeBlock>());

    params.setBlockValue("GAUGE", 2, std::pow(2, 1.25) * m_W * std::sqrt(G_F));  // g2
    params.setBlockValue("GAUGE", 1, params("GAUGE", 2) * sW / std::sqrt(1 - sW * sW)); // gp 
    params.setBlockValue("GAUGE", 3, std::sqrt(4 * M_PI * alpha_s_MZ)); // gs_MZ
    params.setBlockValue("GAUGE", 4, std::sqrt(4 * M_PI / inv_alpha_em)); // e_em

    std::string assets_path = project_assets_root.data();;
    JSONParser::getInstance(0)->saveToFile(assets_path + "savestate/parameters_SM.json");
}

// SUSYModelStrategy implementation
void SUSYModelStrategy::initializeParameters(Parameters& params) {

    auto mm = MemoryManager::GetInstance();
    LhaReader* lha = mm->getReader();

    if (!mm->isSpectrum()) {
        std::string root = project_assets_root.data();
        std::string spectrumFile = root + "spectrum/" + mm->getInputLhaPath().filename().c_str();
        LOG_INFO("Starting SUSY spectrum calculation...");
        CalculatorType calculatorType = CalculatorType::Softsusy;
        GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", mm->getInputLhaPath(), spectrumFile);
        LOG_INFO("THDM spectrum calculation ran sucessfully");
        lha->update(spectrumFile);
        LOG_INFO("LHA Blocks updated.");
    }

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
        std::shared_ptr<ArrayBlock<2, 2>> uniquePtr = std::make_shared<ArrayBlock<2,2>>(std::move(*elem.second));
        uniquePtr->setValues(temp);
        params.addBlock(elem.first, uniquePtr);
    }

    for (auto& elem : type44) {
        std::array<std::array<double, 4>,4> temp;
        readMatrix(temp, elem.first, lha);
        std::shared_ptr<ArrayBlock<4, 4>> uniquePtr = std::make_shared<ArrayBlock<4,4>>(std::move(*elem.second));
        uniquePtr->setValues(temp);
        params.addBlock(elem.first, uniquePtr);
    }

    for (auto& elem : type33) {
        std::array<std::array<double, 3>,3> temp;
        readMatrix(temp, elem.first, lha);
        std::shared_ptr<ArrayBlock<3, 3>> uniquePtr = std::make_shared<ArrayBlock<3,3>>(std::move(*elem.second));
        uniquePtr->setValues(temp);
        params.addBlock(elem.first, uniquePtr);
    }

    auto alphablock = std::make_shared<AlphaBlock>();
    alphablock->setValue(0, lha->getValue<double>("ALPHA", 0));
    params.addBlock("ALPHA", std::move(alphablock));

    auto hmixblock = std::make_shared<HMIXBlock>();
    hmixblock->setValue(0, static_cast<LhaElement<double>*>(lha->getBlock("HMIX")->get(1))->getScale());
    std::vector<double> values (4);
    lha->extractFromBlock("HMIX", values);
    for (size_t i=0; i!=values.size(); ++i) {
        hmixblock->setValue(i+1, values[i]);
    }
    params.addBlock("HMIX", std::move(hmixblock));

    values.resize(3);
    auto gaugeblock = std::make_shared<GaugeBlock>();
    lha->extractFromBlock("GAUGE", values);
    for (size_t i=0; i!=values.size(); ++i) {
        gaugeblock->setValue(i+1, values[i]);
    }
    params.addBlock("GAUGE", std::move(gaugeblock));

    auto msoftblock = std::make_shared<MSOFTBlock>();
    auto elts = lha->getBlock("MSOFT")->getEntries();
    for (size_t i = 0; i < elts->size(); ++i) {
        auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
        msoftblock->setValue(e->getId(), e->getValue());
    }
    params.addBlock("MSOFT", std::move(msoftblock));

    auto massblock = std::make_shared<MassBlock>();
    elts = lha->getBlock("MASS")->getEntries();
    massblock->setValue(1000039, 0.);
    massblock->setValue(45, 0.);
    massblock->setValue(46, 0.);
    for (size_t i = 0; i < elts->size(); ++i) {
        auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
        massblock->setValue(e->getId(), e->getValue());
    }
    params.addBlock("MASS", std::move(massblock));

    std::string root_path = project_assets_root.data();
    JSONParser::getInstance(0)->saveToFile(root_path+ "savestate/parameters_SUSY.json");
}

// THDMModelStrategy implementation
void THDMModelStrategy::initializeParameters(Parameters& params) {
    LhaReader* lha = MemoryManager::GetInstance()->getReader();
    MemoryManager * memo = MemoryManager::GetInstance();

    if (!memo->isSpectrum()) {
        std::string root = project_assets_root.data();
        std::string spectrumFile = root + "spectrum/" + memo->getInputLhaPath().filename().c_str();
        LOG_INFO("Starting THDM spectrum calculation...");
        CalculatorType calculatorType = CalculatorType::TwoHDM;
        GeneralCalculatorFactory::executeCommand(calculatorType, "calculateSpectrum", memo->getInputLhaPath(), spectrumFile);
        LOG_INFO("THDM spectrum calculation ran sucessfully");
        lha->update(spectrumFile);
        LOG_INFO("LHA Blocks updated.");
    }

    std::vector<std::string> mandatory {"MINPAR", "MASS", "ALPHA"};

    auto massblock = std::make_shared<MassBlock>();
    auto elts = lha->getBlock("MASS")->getEntries();
    for (size_t i = 0; i < elts->size(); ++i) {
        auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
        massblock->setValue(e->getId(), e->getValue());
        // this->masses[std::stoi(e->getId())] = e->getValue();
    }
    params.addBlock("MASS", std::move(massblock));

    auto alphablock = std::make_shared<AlphaBlock>();
    alphablock->setValue(0, lha->getValue<double>("ALPHA", 0));
    params.addBlock("ALPHA", std::move(alphablock));
    auto hmixblock = std::make_shared<HMIXBlock>();
    double tan_beta = lha->getValue<double>("MINPAR", 3);
    hmixblock->setValue(2, tan_beta);
    params.addBlock("HMIX", std::move(hmixblock));
    // std::cout << "dd : " << lha->getValue<int>("MINPAR", "24") << std::endl;
     //small patch
     int type;
     try {
        type = static_cast<int>(lha->getValue<double>("MINPAR", 24));
     } catch(...) {
        type = 4;
     }
    double cot_beta = 1 / tan_beta;

    auto yublock = std::make_shared<YUBlock>();
    auto ydblock = std::make_shared<YUBlock>();
    auto yeblock = std::make_shared<YUBlock>();

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
            yeblock->setValue(11, -tan_beta);
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

    std::string root_path = project_assets_root.data();
    JSONParser::getInstance(0)->saveToFile(root_path+ "savestate/parameters_THDM.json");

}

// FlAVORModelStrategy implementation
void FlavorStrategy::initializeParameters(Parameters& params) {
    auto mm = MemoryManager::GetInstance();
    LhaReader* lha = mm->getReader();
    std::vector<std::string> blocks = lha->getBlocksNames();
    std::vector<Correlation> _;
    std::vector<Value> values;
    read_json(mm->getParameterCovariancePath().string(), values, _);

    auto massblock = std::make_shared<FMassBlock>();
    populate_from_json(massblock, values);
    overwrite_from_lha(massblock, lha);
    params.addBlock("FMASS", std::move(massblock));

    auto lifetimeblock = std::make_shared<FLifeBlock>();
    populate_from_json(lifetimeblock, values);
    overwrite_from_lha(lifetimeblock, lha);
    params.addBlock("FLIFE", std::move(lifetimeblock));

    auto fconstblock = std::make_shared<FConstBlock>();
    for (const auto& val : values) {
        if (val.name.starts_with("FCONST")) {
            std::string del = "|";
            auto split = val.name.find(del);
            std::string pdg_str = val.name.substr(split + 1, val.name.size() - split);
            split = pdg_str.find(del);
            int pdg = 100 * std::stoi(pdg_str.substr(0, split)) + std::stoi(pdg_str.substr(split + 1, pdg_str.size() - split));
            fconstblock->setValue(pdg, val.central_value);
        }
    }
    if (std::find(blocks.begin(), blocks.end(), "FCONST") != blocks.end()) {
        fconstblock->setValue(51101, lha->getValue<double>("FCONST", {511, 1})); // f_B
        fconstblock->setValue(52101, lha->getValue<double>("FCONST", {521, 1})); // f_B0
        fconstblock->setValue(53101, lha->getValue<double>("FCONST", {531, 1})); // f_Bs
        fconstblock->setValue(32301, lha->getValue<double>("FCONST", {323, 1})); // f_K*_par
        fconstblock->setValue(32302, lha->getValue<double>("FCONST", {323, 2})); // f_K*_perp
    }
    params.addBlock("FCONST", std::move(fconstblock));
}   

void GeneralModelStrategy::initializeParameters(Parameters& params) {
    LhaReader* lha = MemoryManager::GetInstance()->getReader();
    MemoryManager * memo = MemoryManager::GetInstance();
    std::string spectrumFile = MemoryManager::GetInstance()->getInputLhaPath();

    for (auto elem : lha->getBlocksNames()) {
        LOG_DEBUG(elem);
        LOG_DEBUG(lha->findPrototype(elem).itemCount);
        if (lha->findPrototype(elem).itemCount == 2) {
            auto generalblock = std::make_shared<MapBlock>();
            generalblock->blockname = elem;
            auto elts = lha->getBlock(elem)->getEntries();
            for (size_t i = 0; i < elts->size(); ++i) {
                auto e = static_cast<LhaElement<double>*>(elts->at(i).get());
                generalblock->setValue(e->getId(), e->getValue());
            }
            params.addBlock(elem, std::move(generalblock));
            continue;
        }
        
    }
    lha->update(spectrumFile);
    
    LOG_INFO("LHA Blocks updated.");

    std::vector<std::string> mandatory {"MINPAR", "MASS"};

    std::string root_path = project_assets_root.data();
    JSONParser::getInstance(0)->saveToFile(root_path+ "savestate/parameters_GENERAL.json");
}

void WilsonInputStrategy::initializeParameters(Parameters &params) {
    auto fill_wilson_block = [] (const std::string& block_name, const std::string& flha_name, double scale, int type, Parameters& params, std::vector<int>& nonzero) -> std::pair<double, int> {
        params.addBlock(block_name, std::make_shared<WilsonBlock>());
        auto lha = MemoryManager::GetInstance()->getReader();
        auto block = lha->getBlock(flha_name);
        if (!block) {
            if (flha_name == "FWCOEF")
                LOG_ERROR("Parameters", "Unable to find real parts of wilson coefficients (block FWCOEF not found in FLHA file)");
            return {scale, type};
        }
        for (auto &e : *(block->getEntries())) {
            int content = e->getId().parts[0];
            int structure = e->getId().parts[1];
            int order = e->getId().parts[2];

            if (order >= 10) {
                LOG_WARN("Found QED corrections to Wilson coefficient, skipping");
                continue;
            }
            
            auto c = static_cast<LhaElement<double>*>(e.get());
            if (scale == -1)
                scale = c->getScale();
            else if (scale != c->getScale()) {
                LOG_ERROR("Parameters", "All Wilson coefficients must be given at the same scale.");
            }

            if (type == -1)
                type = e->getId().parts[3];
            else if (type != e->getId().parts[3]) {
                LOG_ERROR("Parameters", "All Wilson coefficients must be of the same type.");
            }

            params.setBlockValue(block_name, 10 * (int)WCoefMapper::from_flha(content, structure) + order, c->getValue());
            nonzero.push_back((int)WCoefMapper::from_flha(content, structure));
        }

        params.setBlockValue(block_name, -1, scale);
        params.setBlockValue(block_name, -2, type);
        return {scale, type};
    };

    if (!MemoryManager::GetInstance()->hasWilsons()) {
        LOG_ERROR("Parameters", "No Wilson coefficients were given in the input file.");
    }

    std::vector<int> nonzero_re = {};
    std::vector<int> nonzero_im = {};
    double scale = -1;
    int type = -1;
    auto p = fill_wilson_block("REWCOEF", "FWCOEF", scale, type, params, nonzero_re);
    scale = p.first;
    type = p.second;
    fill_wilson_block("IMWCOEF", "IMFWCOEF", scale, type, params, nonzero_im);

    for (int i = 0; i < WCoefMapper::n_wilsons(); i++) {
        if (std::find(nonzero_re.begin(), nonzero_re.end(), i) == nonzero_re.end()) {
            params.setBlockValue("REWCOEF", i * 10, 0);
        }
        if (std::find(nonzero_im.begin(), nonzero_im.end(), i) == nonzero_im.end()) {
            params.setBlockValue("IMWCOEF", i * 10, 0);
        }
    }   
}

void FormFactorStrategy::initializeParameters(Parameters &params) {
    std::vector<Correlation> _;
    std::vector<Value> values;
    read_json(MemoryManager::GetInstance()->getParameterCovariancePath().string(), values, _);

    std::vector<std::shared_ptr<MapBlock>> ff_blocks = {
        std::make_shared<BKsBlock>(),
        std::make_shared<BllBlock>(),
        std::make_shared<BXsBlock>(),
        std::make_shared<BDlnuBlock>(),
        std::make_shared<BDslnuBlock>(),
    };

    for (auto b : ff_blocks) {
        populate_from_json(b, values);
        params.addBlock(b->blockname, b);
    }
}

void Parameters::changeParameterMode(const ParamId &param_id, ParameterMode new_mode) {
    blockAccessor.setMode(param_id.block, param_id.code, new_mode);
}

void Parameters::shiftParameter(const ParamId &param_id, double shift_value) {
    blockAccessor.setValue(param_id.block, param_id.code, blockAccessor.getValue(param_id.block, param_id.code) + shift_value, true);
}

std::shared_ptr<Parameters> ParametersFactory::GetParameters(ParameterType id) {
    if (instances.find(id) == instances.end()) {
        std::shared_ptr<ModelStrategy> strategy = createStrategy(id);
        instances[id] = std::make_shared<Parameters>(Parameters(strategy));
    }
    return instances[id];
}

void ParametersFactory::removeParameters(ParameterType id) {
    if (instances.find(id) == instances.end()) {
        LOG_ERROR("OutOfRange", "Cannot remove parameters if it doesn't exist");
    }
    LOG_DEBUG("erasing ; ", (int)id);
    std::shared_ptr<Parameters> _ = instances[id];
    instances.erase(id);
}

std::shared_ptr<ModelStrategy> ParametersFactory::createStrategy(ParameterType id) {
    switch (id) {
        case ParameterType::SM:
            return std::make_shared<SMModelStrategy>(SMModelStrategy());
        case ParameterType::SUSY:
            return std::make_shared<SUSYModelStrategy>(SUSYModelStrategy());
        case ParameterType::THDM:
            return std::make_shared<THDMModelStrategy>(THDMModelStrategy());
        case ParameterType::FLAVOR:
            return std::make_shared<FlavorStrategy>(FlavorStrategy());
        case ParameterType::CUSTOM:
            return std::make_shared<GeneralModelStrategy>(GeneralModelStrategy());
        case ParameterType::WILSON:
            return std::make_shared<WilsonInputStrategy>(WilsonInputStrategy());
        case ParameterType::FF:
            return std::make_shared<FormFactorStrategy>(FormFactorStrategy());
        default:
            throw std::invalid_argument("Unknown parameters instance ID");
    }
}
