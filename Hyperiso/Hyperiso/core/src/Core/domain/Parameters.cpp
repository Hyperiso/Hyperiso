#include "Parameters.h"
#include "DependentBlockManager.h"

std::map<ParameterType, std::shared_ptr<Parameters>> Parameters::instances;
std::map<ParameterType, std::shared_ptr<Parameters>> ParametersFactory::instances;

std::shared_ptr<Parameters> Parameters::GetInstance(ParameterType id) {
    LOG_DEBUG("Trying to access Parameters instance of type", static_cast<int>(id));
    auto allowed = MemoryManager::GetInstance()->getMemoryCache().parameter_types;
    if (std::find(allowed.begin(), allowed.end(), id) == allowed.end())
        LOG_ERROR("OutOfRange", "Parameter type undefined: ", ParameterTypeMapper::str(id));
    return ParametersFactory::GetParameters(id);
}

void Parameters::CleanupInstance(ParameterType id) {
    ParametersFactory::removeParameters(id);
}

void Parameters::claim_parameters(ParameterType type) {
    for (auto& [_, block] : *blockAccessor) {
        block->set_owner(type);
    }
}

Parameters::Parameters(std::shared_ptr<ModelStrategy> modelStrategy)
    : strategy(modelStrategy)
{
    LOG_VERBOSE("Param creation at", this);
    auto truc = strategy->initializeParameters(*this);
    strategy->add_absent_block(truc);
}

scalar_t Parameters::operator()(const std::string& block, LhaID id) const {
    // if (block == "WPARAM_MATCH_SM" ){

    //     std::cout << blockAccessor << std::endl;
    // }
    return blockAccessor->getValue(block, id);
}

std::shared_ptr<Parameter> Parameters::get_parameter(const std::string &block,
                                                     LhaID pdgCode) {
    return blockAccessor->at(block)->retrieve(pdgCode);
}

bool Parameters::exist(const std::string& block, LhaID id) {
    return blockAccessor->has_param(block, id);
}

void Parameters::setBlockValue(const std::string& name, LhaID id, double value) {
    blockAccessor->setValue(name, id, value);
}

std::map<LhaID, double> Parameters::get_block_infos(std::string blockName) {
    return blockAccessor->getAllValues(blockName);
}

std::unordered_set<std::string> Parameters::get_blocks_list() {
    return blockAccessor->get_block_names();
}

std::unordered_set<std::string> Parameters::init_blocks(ParameterType type) {
    std::unordered_set<std::string> existing, missing;
    
    std::ranges::partition_copy(
        ParameterBlockRepartition::BLOCKS.at(type),
        std::inserter(existing, existing.end()),
        std::inserter(missing, missing.end()),
        [&](const std::string& s) {
            return MemoryManager::GetInstance()->input_cache->get_block_names().contains(s);
        }
    );

    if (type == ParameterType::WILSON && !MemoryManager::GetInstance()->cache.config.flags[ExternalFlag::HAS_WILSON_INPUT]) {
        existing.erase("FWCOEF");
        existing.erase("IMFWCOEF");
    }

    if (type == ParameterType::OBSERVABLE && !MemoryManager::GetInstance()->cache.config.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT]) {
        return missing;
    }
    
    this->blockAccessor = MemoryManager::GetInstance()->extract_blocks(existing);
    claim_parameters(type);
    return missing;
}

std::unordered_set<std::string> SMModelStrategy::initializeParameters(Parameters& params) {
    
    auto absent_blocks = params.init_blocks(ParameterType::SM);


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
    std::cout << "fuuuuck strategy" << std::endl;
    

    // TODO : Initialize derived blocks RE/IMUPMNS
    // TODO : Calculate W mass and store it somewhere
    // TODO : Export savestate to JSON
    return absent_blocks;
}

void SMModelStrategy::postInitialization(Parameters& params) {
    QCDHelper::Init();

    if (absent_blocks.contains("VCKM")) {
        std::unordered_map<ParameterType, std::vector<std::string>> src = {
            {ParameterType::SM, {"VCKMIN"}}
        };
    
        auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
            std::cout << "ITS MEEEE, MARIOOO" << std::endl;
            double lambda = src.at("VCKMIN")->retrieve(1)->get_val();
            double l2 = lambda * lambda;
            double l3 = l2 * lambda;
            double A = src.at("VCKMIN")->retrieve(2)->get_val();
            double rho = src.at("VCKMIN")->retrieve(3)->get_val();
            double eta = src.at("VCKMIN")->retrieve(4)->get_val();

            dep_block->store_or_assign(LhaID(0, 0), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(0, 0)}, 1 - l2 / 2, 0., 0.));
            dep_block->store_or_assign(LhaID(0, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(0, 1)}, lambda, 0., 0.));
            dep_block->store_or_assign(LhaID(0, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(0, 2)}, A * l3 * complex_t(rho, -eta), 0., 0.));
            dep_block->store_or_assign(LhaID(1, 0), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(1, 0)}, -lambda, 0., 0.));
            dep_block->store_or_assign(LhaID(1, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(1, 1)}, 1 - l2 / 2, 0., 0.));
            dep_block->store_or_assign(LhaID(1, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(1, 2)}, A * l2, 0., 0.));
            dep_block->store_or_assign(LhaID(2, 0), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(2, 0)}, A * l3 * complex_t(1 - rho, -eta), 0., 0.));
            dep_block->store_or_assign(LhaID(2, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(2, 1)}, -A * l2, 0., 0.));
            dep_block->store_or_assign(LhaID(2, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(2, 2)}, 1, 0., 0.));
        };

        DependentBlockManager::addDependentBlock("VCKM", src, ParameterType::SM, func);
    }
}

std::unordered_set<std::string> BSMModelStrategy::initializeParameters(Parameters& params) {
    auto absent_blocks = params.init_blocks(ParameterType::BSM);
    return absent_blocks;
    // TODO : Export savestate to JSON
}

// void SUSYModelStrategy::initializeParameters(Parameters& params) {
//     params.init_blocks(ParameterType::BSM);
//     // TODO : Export savestate to JSON
// }

// void THDMModelStrategy::initializeParameters(Parameters& params) {
//     params.init_blocks(ParameterType::BSM);
//     // TODO : Export savestate to JSON
// }

std::unordered_set<std::string> FlavorStrategy::initializeParameters(Parameters& params) {
    auto absent_blocks = params.init_blocks(ParameterType::FLAVOR);
    return absent_blocks;
    // TODO : Export savestate to JSON
}   

// void GeneralModelStrategy::initializeParameters(Parameters& params) {
//     params.init_blocks(ParameterType::BSM);
//     // TODO : Export savestate to JSON
// }

std::unordered_set<std::string> WilsonInputStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::WILSON);
    return absent_blocks;
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

std::unordered_set<std::string> DecayStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::DECAY);
    return absent_blocks;
    // TODO : Export savestate to JSON
}

std::unordered_set<std::string> ObservableStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::OBSERVABLE);
    return absent_blocks;
}

std::unordered_set<std::string> PassthroughStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::PASSTHROUGH);
    return absent_blocks;
}


void Parameters::changeParameterMode(const ParamId &param_id, ParameterMode new_mode) {
    // blockAccessor->setMode(param_id.block, param_id.code, new_mode);
}

void Parameters::shiftParameter(const ParamId &param_id, double shift_value) {
    blockAccessor->setValue(param_id.block, param_id.code, blockAccessor->getValue(param_id.block, param_id.code) + shift_value);
}

std::shared_ptr<Parameters> ParametersFactory::GetParameters(ParameterType id) {
    if (instances.find(id) == instances.end()) {
        std::shared_ptr<ModelStrategy> strategy = createStrategy(id);
        instances[id] = std::make_shared<Parameters>(Parameters(strategy));
        strategy->postInitialization(*instances[id]);
        std::cout << (int)id << std::endl;
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
        case ParameterType::BSM:
            return std::make_shared<BSMModelStrategy>(BSMModelStrategy());
        case ParameterType::FLAVOR:
            return std::make_shared<FlavorStrategy>(FlavorStrategy());
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