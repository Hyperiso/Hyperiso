/**
 * @file Parameters.h
 * @brief Defines strategies for different physics models and parameter management.
 * 
 * This file declares several Parameters singletons instances using strategy classes that manage initialization and
 * manipulation of the instances.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "BlockAccessor.h"
#include "MemoryManager.h"
#include "Interface.h"
#include "QCDHelper.h"
#include "config.hpp"
#include <memory>
#include <ranges>
#include <algorithm>
#include <unordered_set>
#include "ParameterRouter.h"

/**
 * @class ModelStrategy
 * @brief Abstract base class for model-specific strategies.
 */
class ModelStrategy {
public:
    /**
     * @brief Initializes model-specific parameters.
     * @param params Reference to Parameters object.
     */
    virtual std::unordered_set<std::string> initializeParameters(class Parameters& params) = 0;
    virtual void postInitialization(Parameters& params) = 0;
    virtual ~ModelStrategy() = default;
    void add_absent_block(std::unordered_set<std::string> _) {absent_blocks = _;};
    void remove_absent_block() {absent_blocks = std::unordered_set<std::string>();}
protected:
    std::unordered_set<std::string> absent_blocks;
};

// Concrete Strategy Classes
/** @class SMModelStrategy @brief Strategy for the Standard Model. */
class SMModelStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override;
};

/** @class BSMModelStrategy @brief Strategy for BSM models. */
class BSMModelStrategy : public ModelStrategy {
    public:
    std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
        void postInitialization(Parameters& params) override {}
    };

// /** @class SUSYModelStrategy @brief Strategy for SUSY models. */
// class SUSYModelStrategy : public ModelStrategy {
// public:
//     void initializeParameters(class Parameters& params) override;
// };

// /** @class THDMModelStrategy @brief Strategy for Two-Higgs-Doublet Models. */
// class THDMModelStrategy : public ModelStrategy {
// public:
//     void initializeParameters(class Parameters& params) override;
// };

/** @class FlavorStrategy @brief Strategy for flavor physics parameters. */
class FlavorStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** @class GeneralModelStrategy @brief General model strategy for parameter initialization. */
class GeneralModelStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** @class WilsonInputStrategy @brief Strategy for Wilson coefficient inputs. */
class WilsonInputStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** @class DecayStrategy @brief Strategy for decay parameters (incl. formfactors). */
class DecayStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** @class ObservableStrategy @brief Strategy for observable inputs (exp. values). */
class ObservableStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

    /** @class PassthroughStrategy @brief Strategy for passthrough parameters (not needed at runtime but must be in the output). */
class PassthroughStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/**
 * @class Parameters
 * @brief Manages parameter values and strategies for different models. Manage every parameters stored
 * which came from the lha.
 */
class Parameters {
public:
    /**
     * @brief Retrieves an instance of Parameters for a given model type.
     * @param id The model type (default: Standard Model).
     * @return Shared pointer to a Parameters instance.
     */
    static std::shared_ptr<Parameters> GetInstance(ParameterType id = ParameterType::SM);

    /**
     * @brief Cleans up and removes an instance of Parameters for a given model type.
     * @param id The parameters type to remove.
     */
    void CleanupInstance(ParameterType id = ParameterType::SM);

    /**
     * @brief Checks if a parameter exists within a specified block.
     * @param block The name of the block.
     * @param pdgCode The PDG code to check.
     * @return True if the parameter exists, false otherwise.
     */
    bool exist(const std::string& block, LhaID pdgCode);
    
    /**
     * @brief Retrieves a parameter value using function call syntax.
     * @param block The name of the block.
     * @param pdgCode The PDG code.
     * @return The corresponding parameter value.
     */
    scalar_t operator()(const std::string& block, LhaID pdgCode) const;

    std::shared_ptr<Parameter> get_parameter(const std::string& block, LhaID pdgCode);

    /**
     * @brief Sets a parameter value within a specified block.
     * @param name The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to assign.
     * @param force If true, forces the update.
     */
    void setBlockValue(const std::string& name, LhaID pdgCode, double value);

    /**
     * @brief Retrieves all parameter values from a specified block.
     * @param blockName The name of the block.
     * @return A map of PDG codes to parameter values.
     */
    std::map<LhaID, double> get_block_infos(std::string blockName);

    /**
     * @brief Retrieves a list of all available parameter blocks.
     * @return A set of block names.
     */
    std::unordered_set<std::string> get_blocks_list();

    /**
     * @brief Changes the operational mode of a specified parameter.
     * @param param_id The unique identifier of the parameter.
     * @param new_mode The new mode to apply.
     */
    void changeParameterMode(const ParamId& param_id, ParameterMode new_mode);

    /**
     * @brief Adjusts the value of a specified parameter by a shift amount.
     * @param param_id The unique identifier of the parameter.
     * @param shift_value The amount by which to shift the parameter value.
     */
    void shiftParameter(const ParamId& param_id, double shift_value);

    std::unordered_set<std::string> init_blocks(ParameterType type);

    void freeze_block(const std::string& blockName);
    void unfreeze_block(const std::string& blockName);

    void freeze_param(const std::string& blockName, const LhaID& id);
    void unfreeze_param(const std::string& blockName, const LhaID& id);

    /**
     * @brief Destructor that logs the destruction of a Parameters instance.
     */
    ~Parameters() { LOG_DEBUG("Parameters at ", this); }

private:
    void claim_parameters(ParameterType type);
    
    /** @brief Private constructor for singleton pattern. */
    explicit Parameters(std::shared_ptr<ModelStrategy> modelStrategy);

    /** @brief Map of model types to their corresponding Parameters instances. */
    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;

    /** @brief Strategy used for parameter management. */
    std::shared_ptr<ModelStrategy> strategy;

    std::shared_ptr<BlockAccessor> blockAccessor;

    /** @brief Factory friend. */
    friend class ParametersFactory;
    friend class DependentBlockManager;
};

/**
 * @class ParametersFactory
 * @brief Factory class for creating and managing Parameters instances.
 * 
 * The factory is responsible for ensuring that parameter instances are created
 * and managed correctly. It provides methods to retrieve and remove instances
 * of Parameters for different models.
 */
class ParametersFactory {
public:
    /**
     * @brief Retrieves an instance of Parameters for a given model type.
     * @param id The type of model parameters to retrieve.
     * @return A shared pointer to the Parameters instance.
     */
    static std::shared_ptr<Parameters> GetParameters(ParameterType id);

    /**
     * @brief Removes an instance of Parameters for a given model type.
     * @param id The model type whose instance should be removed.
     */
    static void removeParameters(ParameterType id);
private:
    /**
     * @brief Stores instances of Parameters for different model types.
     */
    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;

    /**
     * @brief Creates an appropriate ModelStrategy based on the given model type.
     * @param id The model type.
     * @return A shared pointer to the newly created ModelStrategy.
     */
    static std::shared_ptr<ModelStrategy> createStrategy(ParameterType id);
};

#endif