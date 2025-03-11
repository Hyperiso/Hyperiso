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

void Parameters::addDependantBlock(const std::string& name, std::shared_ptr<DependentBlock>& block, const std::string& source_block, std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateFunc) {
    blockAccessor->addDependentBlock(name, block, source_block, recalculateFunc);
}

void Parameters::setBlockValue(const std::string& name, LhaID id, double value, bool force) {
    blockAccessor->setValue(name, id, value, force);
}

std::map<LhaID, double> Parameters::get_block_infos(std::string blockName) {
    return blockAccessor->getAllValues(blockName);
}

std::unordered_set<std::string> Parameters::get_blocks_list() {
    return blockAccessor->get_block_names();
}

std::unordered_set<std::string> Parameters::init_blocks(ParameterType type) {
    if (type == ParameterType::CUSTOM) {
        auto block_names = MemoryManager::GetInstance()->get_all_blocks();
        this->blockAccessor = MemoryManager::GetInstance()->extract_blocks(ParamRouter::GetOwnedBlocks(type));
    } else {
        std::unordered_set<std::string> existing, missing;
        std::ranges::partition_copy(
            ParameterBlockRepartition::BLOCKS.at(type),
            std::inserter(existing, existing.end()),
            std::inserter(missing, missing.end()),
            [&](const std::string& s) {
                return MemoryManager::GetInstance()->get_all_blocks().contains(s);
            }
        );
        
        this->blockAccessor = MemoryManager::GetInstance()->extract_blocks(existing);
        return missing;
    }
}

void SMModelStrategy::initializeParameters(Parameters& params) {
    
    auto absent_blocks = params.init_blocks(ParameterType::SM);
    QCDHelper::Init(params("SMINPUTS", 3), params("SMINPUTS", 4), params("SMINPUTS", 6), params("SMINPUTS", 5),  
                    params("MASS", 4), params("MASS", 3), params("MASS", 2), params("MASS", 1));

    std::shared_ptr<DependentBlock> gauge_block = nullptr;

    auto gauge_update_func = [](std::shared_ptr<Block> src, std::shared_ptr<DependentBlock> dep_block) {
        double e_em = std::sqrt(4 * PI / src->getValue(1));
        double g_3 = std::sqrt(4 * PI * src->getValue(3));
        // Add loop corrections to theta_w
        double theta_w = 0.5 * std::asin(std::sqrt(4 * PI * INV_RT2 / (src->getValue(1) * src->getValue(2))) / src->getValue(4));
        dep_block->setValue(1, e_em / std::sin(theta_w), true);
        dep_block->setValue(2, e_em / std::cos(theta_w), true);
        dep_block->setValue(3, std::sqrt(4 * PI * src->getValue(3)));
        dep_block->setValue(4, e_em);
    };

    params.addDependantBlock("GAUGE", gauge_block, "SMINPUTS", gauge_update_func);

    // std::shared_ptr<DependentBlock> gauge_block = nullptr;

    // auto gauge_update_func = [](std::shared_ptr<Block> src, std::shared_ptr<DependentBlock> dep_block) {
    //     double e_em = std::sqrt(4 * PI / src->getValue(1));
    //     double g_3 = std::sqrt(4 * PI * src->getValue(3));
    //     // Add loop corrections to theta_w
    //     double theta_w = 0.5 * std::asin(std::sqrt(4 * PI * INV_RT2 / (src->getValue(1) * src->getValue(2))) / src->getValue(4));
    //     dep_block->setValue(1, e_em / std::sin(theta_w), true);
    //     dep_block->setValue(2, e_em / std::cos(theta_w), true);
    //     dep_block->setValue(3, std::sqrt(4 * PI * src->getValue(3)));
    //     dep_block->setValue(4, e_em);
    // };

    // params.addDependantBlock("GAUGE", gauge_block, "SMINPUTS", gauge_update_func);

    // if (absent_blocks.contains("RECKM")) {
    //     std::shared_ptr<ReCKMBlock> reckm_block = nullptr;
    //     params.addDependantBlock("RECKM", reckm_block, "SMINPUTS");
    // }

    // if (absent_blocks.contains("IMCKM")) {
    //     std::shared_ptr<ImCKMBlock> imckm_block = nullptr;
    //     params.addDependantBlock("IMCKM", imckm_block, "SMINPUTS");
    // }

    // TODO : Initialize derived blocks RE/IMUPMNS
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

void DecayStrategy::initializeParameters(Parameters &params) {
    params.init_blocks(ParameterType::DECAY);

    // TODO : Export savestate to JSON
}

void ObservableStrategy::initializeParameters(Parameters &params) {
    params.init_blocks(ParameterType::OBSERVABLE);
}

void PassthroughStrategy::initializeParameters(Parameters &params) {
    params.init_blocks(ParameterType::PASSTHROUGH);
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
            return std::make_shared<DecayStrategy>(DecayStrategy());
        case ParameterType::OBSERVABLE:
            return std::make_shared<DecayStrategy>(DecayStrategy());
        case ParameterType::PASSTHROUGH:
            return std::make_shared<DecayStrategy>(DecayStrategy());
        default:
            throw std::invalid_argument("Unknown parameters instance ID");
    }
}