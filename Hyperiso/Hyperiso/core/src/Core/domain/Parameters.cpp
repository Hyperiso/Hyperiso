#include "Parameters.h"
#include "DependentBlockManager.h"

std::map<ParameterType, std::shared_ptr<Parameters>> Parameters::instances;
std::map<ParameterType, std::shared_ptr<Parameters>> ParametersFactory::instances;

std::shared_ptr<Parameters> Parameters::GetInstance(ParameterType id) {
    // LOG_INFO("Trying to access Parameters instance of type", static_cast<int>(id));
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

scalar_t Parameters::operator()(const BlockName& block, LhaID id) const {
    // if (block == "WPARAM_MATCH_SM" ){

    //     std::cout << blockAccessor << std::endl; 
    // }
    return blockAccessor->getValue(block, id);
}

std::shared_ptr<Parameter> Parameters::get_parameter(const BlockName &block,
                                                     LhaID pdgCode) {
    return blockAccessor->at(block)->retrieve(pdgCode);
}

bool Parameters::exist(const BlockName& block, LhaID id) {
    return blockAccessor->has_param(block, id);
}

void Parameters::setBlockValue(const BlockName& name, LhaID id, scalar_t value) {
    blockAccessor->setValue(name, id, value);
}

std::map<LhaID, scalar_t> Parameters::get_block_infos(BlockName blockName) {
    return blockAccessor->getAllValues(blockName);
}

double Parameters::get_block_scale(BlockName blockName) const{
    return this->blockAccessor->at(blockName)->get_scale();
}

std::unordered_set<BlockName> Parameters::get_blocks_list()
{
    return blockAccessor->get_block_names();
}

std::unordered_set<BlockName> Parameters::init_blocks(ParameterType type) {
    LOG_DEBUG("Init blocks for parameter type", ParameterTypeMapper::str(type));
    std::unordered_set<BlockName> existing, missing;
    
    std::ranges::partition_copy(
        ParameterBlockRepartition::BLOCKS.at(type),
        std::inserter(existing, existing.end()),
        std::inserter(missing, missing.end()),
        [&](const BlockName& s) {
            auto data_ = MemoryManager::GetInstance()->input_cache->get_block_names();
            return std::any_of(
                data_.begin(),
                data_.end(),
                [&](const auto& bn) { return bn == s; }
            );
        }
    );

    if (type == ParameterType::WILSON && !MemoryManager::GetInstance()->cache.config.flags[ExternalFlag::HAS_WILSON_INPUT]) {
        existing.erase("FWCOEF");
        existing.erase("IMFWCOEF");
    }

    // TODO v1.1 : manage case of theoretical obs in lha. For now, assume only exp is given. 
    if (type == ParameterType::OBSERVABLE && MemoryManager::GetInstance()->cache.config.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT]) {
        return missing;
    }
    
    this->blockAccessor = MemoryManager::GetInstance()->extract_blocks(existing);
    // LOG_INFO(this->blockAccessor);
    claim_parameters(type);
    return missing;
}

void Parameters::freeze_block(const BlockName &blockName) {
    if (!blockAccessor->contains(blockName)) {
        // std::cout << blockAccessor << std::endl;
        LOG_INFO(blockAccessor);
        LOG_ERROR("Cannot freeze non-existing block", blockName);
    }

    return this->blockAccessor->at(blockName)->freeze();
}

void Parameters::unfreeze_block(const BlockName &blockName) {
    return this->blockAccessor->at(blockName)->unfreeze();
}

void Parameters::freeze_param(const BlockName &blockName, const LhaID &id) {
    if (!blockAccessor->contains(blockName)) {
        LOG_ERROR("Cannot freeze dependent parameter in non-existing block", blockName);
    }

    if (!blockAccessor->at(blockName)->contains(id)) {
        LOG_ERROR("Cannot freeze non-existing dependent parameter", id.to_string(), "in block", blockName);
    }

    return this->blockAccessor->at(blockName)->retrieve(id)->freeze();
}

void Parameters::unfreeze_param(const BlockName &blockName, const LhaID &id) {
    return this->blockAccessor->at(blockName)->retrieve(id)->unfreeze();
}

std::ostream &operator<<(std::ostream &os, std::shared_ptr<Parameters> instance) {
    os << instance->blockAccessor;
    return os;
}

std::unordered_set<BlockName> SMModelStrategy::initializeParameters(Parameters& params) {
    auto absent_blocks = params.init_blocks(ParameterType::SM);
    return absent_blocks;
}

void SMModelStrategy::postInitialization(Parameters& params) {
    // Ask Nazila : Which value to use ? BSM (i.e. from spectrum ?) or SMINPUTS-derived (dependent block) ?
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
    

    // TODO : Initialize derived blocks RE/IMUPMNS
    // Ask Nazila : Calculate W mass and store it somewhere

    QCDHelper::Init();

    if (absent_blocks.contains("VCKM")) {
        std::unordered_map<ParameterType, std::vector<std::string>> src_ckm = {
            {ParameterType::SM, {"VCKMIN"}}
        };
    
        auto func_ckm = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
            double lambda = src.get_val("VCKMIN", 1);
            double l2 = lambda * lambda;
            double l3 = l2 * lambda;
            double A = src.get_val("VCKMIN", 2);
            double rho = src.get_val("VCKMIN", 3);
            double eta = src.get_val("VCKMIN", 4);

            double s_12 = lambda;
            double s_23 = A * l2;
            complex_t u_13 = A * l3 * complex_t(rho, eta) * sqrt(1 - std::pow(A * l2, 2)) / (std::sqrt(1 - l2) * (1. - std::pow(A * l2, 2) * complex_t(rho, eta)));
            double s_13 = std::abs(u_13);
            complex_t expid = u_13 / std::abs(u_13);
            double c_12 = std::sqrt(1 - s_12 * s_12);
            double c_23 = std::sqrt(1 - s_23 * s_23);
            double c_13 = std::sqrt(1 - s_13 * s_13);

            dep_block->store_or_assign(LhaID(0, 0), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(0, 0)}, c_12 * c_13, 0., 0.));
            dep_block->store_or_assign(LhaID(0, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(0, 1)}, s_12 * c_13, 0., 0.));
            dep_block->store_or_assign(LhaID(0, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(0, 2)}, s_13 / expid, 0., 0.));
            dep_block->store_or_assign(LhaID(1, 0), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(1, 0)}, -s_12 * c_23 - c_12 * s_23 * s_13 * expid, 0., 0.));
            dep_block->store_or_assign(LhaID(1, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(1, 1)}, c_12 * c_23 - s_12 * s_23 * s_13 * expid, 0., 0.));
            dep_block->store_or_assign(LhaID(1, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(1, 2)}, s_23 * c_13, 0., 0.));
            dep_block->store_or_assign(LhaID(2, 0), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(2, 0)}, s_12 * s_23 - c_12 * c_23 * s_13 * expid, 0., 0.));
            dep_block->store_or_assign(LhaID(2, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(2, 1)}, -c_12 * s_23 - s_12 * c_23 * s_13 * expid, 0., 0.));
            dep_block->store_or_assign(LhaID(2, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "VCKM", LhaID(2, 2)}, c_23 * c_13, 0., 0.));
        };

        DependentBlockManager::addDependentBlock("VCKM", src_ckm, ParameterType::SM, func_ckm);
    }

    Parameters::GetInstance(ParameterType::WILSON);
    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::SM, {"QCD"}},
        {ParameterType::WILSON, {"EW_SCALE"}}
    };

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        double mu_W = src.get_val("EW_SCALE", 1);
        double mass_top_muW = QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR);
		double mass_b_muW_mbrun = QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR);
		double mass_b_muW_mbpole = QCDHelper::msbar_mass(5, mu_W, MassType::POLE);
		double mass_c_muW = QCDHelper::msbar_mass(4, mu_W, MassType::POLE);

        dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::SM, "MASS_EW_SCALE", 4}, mass_c_muW, 0., 0.));
		dep_block->store_or_assign(LhaID(5, 1), std::make_shared<Parameter>(ParamId{ParameterType::SM, "MASS_EW_SCALE", LhaID(5, 1)}, mass_b_muW_mbrun, 0., 0.));
		dep_block->store_or_assign(LhaID(5, 2), std::make_shared<Parameter>(ParamId{ParameterType::SM, "MASS_EW_SCALE", LhaID(5, 2)}, mass_b_muW_mbpole, 0., 0.));
		dep_block->store_or_assign(6, std::make_shared<Parameter>(ParamId{ParameterType::SM, "MASS_EW_SCALE", 6}, mass_top_muW, 0., 0.));
    };

    DependentBlockManager::addDependentBlock("MASS_EW_SCALE", src, ParameterType::SM, func);
}

std::unordered_set<BlockName> BSMModelStrategy::initializeParameters(Parameters& params) {
    auto absent_blocks = params.init_blocks(ParameterType::BSM);
    return absent_blocks;
}

std::unordered_set<BlockName> FlavorStrategy::initializeParameters(Parameters& params) {
    auto absent_blocks = params.init_blocks(ParameterType::FLAVOR);
    return absent_blocks;
}   

std::unordered_set<BlockName> WilsonInputStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::WILSON);
    return absent_blocks;

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

void WilsonInputStrategy::postInitialization(Parameters& params) {
    if (MemoryManager::GetInstance()->getMemoryCache().config.flags.at(ExternalFlag::HAS_WILSON_INPUT)) {
        params.setBlockValue("EW_SCALE", 1, params.get_block_scale("FWCOEF"));
    }
}

std::unordered_set<BlockName> DecayStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::DECAY);
    return absent_blocks;
    // TODO : Export savestate to JSON
}

std::unordered_set<BlockName> ObservableStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::OBSERVABLE);
    return absent_blocks;
}

std::unordered_set<BlockName> PassthroughStrategy::initializeParameters(Parameters &params) {
    auto absent_blocks = params.init_blocks(ParameterType::PASSTHROUGH);
    return absent_blocks;
}


void Parameters::changeParameterMode(const ParamId &param_id, ParameterMode new_mode) {
    // blockAccessor->setMode(param_id.block, param_id.code, new_mode);
}

void Parameters::shiftParameter(const ParamId &param_id, scalar_t shift_value) {
    blockAccessor->setValue(param_id.block, param_id.code, blockAccessor->getValue(param_id.block, param_id.code) + shift_value);
}

std::shared_ptr<Parameters> ParametersFactory::GetParameters(ParameterType id) {
    if (instances.find(id) == instances.end()) {
        std::shared_ptr<ModelStrategy> strategy = createStrategy(id);
        instances[id] = std::make_shared<Parameters>(Parameters(strategy));
        strategy->postInitialization(*instances[id]);
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
            std::cout << "fck" << std::endl;
            return std::make_shared<ObservableStrategy>(ObservableStrategy());
        case ParameterType::PASSTHROUGH:
            return std::make_shared<PassthroughStrategy>(PassthroughStrategy());
        default:
            throw std::invalid_argument("Unknown parameters instance ID");
    }
}

void Parameters::print_block(const std::string blockname) {
    std::cout << this->blockAccessor->at(blockname) << std::endl;
}