#include "Parameters.h"

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

ParameterType Parameters::GetType(const std::string &block, LhaID pdgCode) {
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

double Parameters::Get(ParameterType type, const std::string& block, LhaID id) {
    LOG_VERBOSE("Attempting to retrieve parameter from instance", (int)type, "with id (", block, ",", id, ")");
    return (*GetInstance(type))(block, id);
}

double Parameters::Get(ParamId id) {
    return Get(id.type, id.block, id.code);
}

Parameters::Parameters(std::shared_ptr<ModelStrategy> modelStrategy)
    : strategy(modelStrategy) { 
    LOG_VERBOSE("Param creation at", this);
    strategy->initializeParameters(*this);
}

double Parameters::operator()(const std::string& block, LhaID id) {
    return blockAccessor->getValue(block, id);
}

bool Parameters::exist(const std::string& block, LhaID id) {
    return blockAccessor->exist(block, id);
}

void Parameters::addBlock(const std::string& name, std::shared_ptr<Block> block) {
    blockAccessor->addBlock(name, block);
}

void Parameters::addDependantBlock(const std::string& name, std::shared_ptr<DependentBlock>& block, const std::string& source_block, std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateFunc) {
    blockAccessor->addDependentBlock(name, block, source_block, recalculateFunc);
}

void Parameters::setBlockValue(const std::string& name, LhaID id, double value, bool force) {
    blockAccessor->setValue(name, id, value, force);
}

std::map<LhaID, double> Parameters::get_block_infos(std::string blockName) {
    return blockAccessor->getAllValues(blockName);
}

std::vector<std::string> Parameters::get_blocks_list() {
    return blockAccessor->get_block_names();
}

complex_t Parameters::get_c_CKM_entry(LhaID id) {
    auto p = Parameters::GetInstance();
    return complex_t((*p)("RECKM", id), (*p)("IMCKM", id));
}

void Parameters::init_blocks(ParameterType type, std::unordered_set<std::string> excluded_dependant) {
    if (type == ParameterType::CUSTOM) {
        auto block_names = MemoryManager::GetInstance()->get_all_blocks();
        this->blockAccessor = MemoryManager::GetInstance()->get_blocks(ParameterBlockRepartition::filter_custom_blocks(block_names));
    } else {
        std::vector<std::string> blocks = ParameterBlockRepartition::BLOCKS.at(type);

        auto it = std::ranges::remove_if(blocks, [&](const std::string& s) {
            return excluded_dependant.contains(s);
        });
        blocks.erase(it.begin(), it.end());
        this->blockAccessor = MemoryManager::GetInstance()->get_blocks(blocks);
    }
}

void SMModelStrategy::initializeParameters(Parameters& params) {
    
    params.init_blocks(ParameterType::SM, {"GAUGE"});
    QCDHelper::Init(params("SMINPUTS", 3), params("SMINPUTS", 4), params("SMINPUTS", 6), params("SMINPUTS", 5),  
                    params("MASS", 4), params("MASS", 3), params("MASS", 2), params("MASS", 1));

    std::cout << "mmh" << std::endl;
    auto recalculateFunc = [](std::shared_ptr<Block> src, std::shared_ptr<DependentBlock> dep_block) {
        std::cout << "Updating fcking block based on " << src->blockname << std::endl;
        double newVal = src->getValue(1) * 2;
        dep_block->setValue(1, newVal, true);
    };

    std::shared_ptr<DependentBlock> gauge_block = nullptr;
    params.addDependantBlock("GAUGE", gauge_block, "SMINPUTS", recalculateFunc);

    params.setBlockValue("SMINPUTS", 1, 10, true);
    // TODO : Initialize derived blocks RE/IMCKM, RE/IMUPMNS
    // TODO : Initialize block GAUGE from SMINPUTS if not already present
    // TODO : Calculate W mass and store it somewhere
    // TODO : Export savestate to JSON
}

void SUSYModelStrategy::initializeParameters(Parameters& params) {
    params.init_blocks(ParameterType::SUSY);

    // TODO : Export savestate to JSON
}

void THDMModelStrategy::initializeParameters(Parameters& params) {
    params.init_blocks(ParameterType::THDM);

    // TODO : Export savestate to JSON
}

void FlavorStrategy::initializeParameters(Parameters& params) {
    params.init_blocks(ParameterType::FLAVOR);

    // TODO : Export savestate to JSON
}   

void GeneralModelStrategy::initializeParameters(Parameters& params) {
    params.init_blocks(ParameterType::CUSTOM);

    // TODO : Export savestate to JSON
}

void WilsonInputStrategy::initializeParameters(Parameters &params) {
    params.init_blocks(ParameterType::CUSTOM);

    // TODO : Export savestate to JSON

    // TODO : Adapt WilsonBlock to MapBlock and rework the following code
    // auto fill_wilson_block = [] (const std::string& block_name, const std::string& flha_name, double scale, int type, Parameters& params, std::vector<int>& nonzero) -> std::pair<double, int> {
    //     params.addBlock(block_name, std::make_shared<WilsonBlock>());
    //     auto lha = MemoryManager::GetInstance()->getReader();
    //     auto block = lha->getBlock(flha_name);
    //     if (!block) {
    //         if (flha_name == "FWCOEF")
    //             LOG_ERROR("Parameters", "Unable to find real parts of wilson coefficients (block FWCOEF not found in FLHA file)");
    //         return {scale, type};
    //     }
    //     for (auto &e : *(block->getEntries())) {
    //         int content = e->getId().parts[0];
    //         int structure = e->getId().parts[1];
    //         int order = e->getId().parts[2];

    //         if (order >= 10) {
    //             LOG_WARN("Found QED corrections to Wilson coefficient, skipping");
    //             continue;
    //         }
            
    //         auto c = static_cast<LhaElement<double>*>(e.get());
    //         if (scale == -1)
    //             scale = c->getScale();
    //         else if (scale != c->getScale()) {
    //             LOG_ERROR("Parameters", "All Wilson coefficients must be given at the same scale.");
    //         }

    //         if (type == -1)
    //             type = e->getId().parts[3];
    //         else if (type != e->getId().parts[3]) {
    //             LOG_ERROR("Parameters", "All Wilson coefficients must be of the same type.");
    //         }

    //         params.setBlockValue(block_name, 10 * (int)WCoefMapper::from_flha(content, structure) + order, c->getValue());
    //         nonzero.push_back((int)WCoefMapper::from_flha(content, structure));
    //     }

    //     params.setBlockValue(block_name, -1, scale);
    //     params.setBlockValue(block_name, -2, type);
    //     return {scale, type};
    // };

    // if (!MemoryManager::GetInstance()->hasWilsons()) {
    //     LOG_ERROR("Parameters", "No Wilson coefficients were given in the input file.");
    // }

    // std::vector<int> nonzero_re = {};
    // std::vector<int> nonzero_im = {};
    // double scale = -1;
    // int type = -1;
    // auto p = fill_wilson_block("REWCOEF", "FWCOEF", scale, type, params, nonzero_re);
    // scale = p.first;
    // type = p.second;
    // fill_wilson_block("IMWCOEF", "IMFWCOEF", scale, type, params, nonzero_im);

    // for (int i = 0; i < WCoefMapper::n_wilsons(); i++) {
    //     if (std::find(nonzero_re.begin(), nonzero_re.end(), i) == nonzero_re.end()) {
    //         params.setBlockValue("REWCOEF", i * 10, 0);
    //     }
    //     if (std::find(nonzero_im.begin(), nonzero_im.end(), i) == nonzero_im.end()) {
    //         params.setBlockValue("IMWCOEF", i * 10, 0);
    //     }
    // }   
}

void FormFactorStrategy::initializeParameters(Parameters &params) {
    params.init_blocks(ParameterType::DECAY);

    // TODO : Export savestate to JSON
}

void Parameters::changeParameterMode(const ParamId &param_id, ParameterMode new_mode) {
    blockAccessor->setMode(param_id.block, param_id.code, new_mode);
}

void Parameters::shiftParameter(const ParamId &param_id, double shift_value) {
    blockAccessor->setValue(param_id.block, param_id.code, blockAccessor->getValue(param_id.block, param_id.code) + shift_value, true);
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
        case ParameterType::DECAY:
            return std::make_shared<FormFactorStrategy>(FormFactorStrategy());
        default:
            throw std::invalid_argument("Unknown parameters instance ID");
    }
}

std::vector<std::string> ParameterBlockRepartition::filter_custom_blocks(const std::vector<std::string> &source) {
    std::vector<std::string> custom_blocks {};

    for (const auto& block : source) {
        if (to_lowercase(block) == "mass" || to_lowercase(block) == "gauge") {
            custom_blocks.emplace_back(block);
            continue;
        }

        for (const auto& [_, known_blocks]: ParameterBlockRepartition::BLOCKS) {
            if (std::find(known_blocks.begin(), known_blocks.end(), block) == known_blocks.end()) {
                custom_blocks.emplace_back(block);
            }
        }
    }

    return custom_blocks;
}
