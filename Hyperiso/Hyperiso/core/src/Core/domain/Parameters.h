/**
 * @file Parameters.h
 * @brief Defines strategies for different physics models and manages parameter instances.
 *
 * This file declares:
 * - ModelStrategy: abstract base class for model-specific strategies.
 * - Concrete strategies (e.g., SMModelStrategy, BSMModelStrategy, etc.).
 * - Parameters: singleton class to manage parameter values and blocks for different models.
 * - ParametersFactory: factory class to create and manage Parameters instances.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <memory>
#include <ranges>
#include <algorithm>
#include <unordered_set>

#include "ParameterRouter.h"
#include "BlockAccessor.h"
#include "MemoryManager.h"
#include "Interface.h"
#include "QCDHelper.h"
#include "config.hpp"

/**
 * @defgroup ParametersModule Parameters and Model Strategies
 * @brief Management of parameters and model-specific initialization.
 *
 * This module defines and manages:
 * - ModelStrategy classes for different physics models (SM, BSM, Flavor, Wilson, etc.)
 * - The Parameters singleton that holds all parameter values.
 * - The ParametersFactory that ensures proper instance creation and destruction.
 *
 * ## Overview
 *
 * The Parameters system relies on a flexible architecture:
 * - Each model type (Standard Model, BSM, etc.) uses a dedicated **ModelStrategy** to load its parameters.
 * - **Parameters** holds the blocks and parameters for each model instance.
 * - **ParametersFactory** guarantees that only one Parameters instance per model exists.
 *
 * ## Structure
 *
 * ```
 * [Input files / MemoryManager]
 *       ↓
 *    [BlockAccessor]
 *       ↓
 * [Parameters instance] ← ModelStrategy
 * ```
 *
 * ## Related Classes
 * - @ref Parameters
 * - @ref ParametersFactory
 * - @ref ModelStrategy
 * - @ref SMModelStrategy
 * - @ref BSMModelStrategy
 * - @ref FlavorStrategy
 * - @ref WilsonInputStrategy
 * - @ref ObservableStrategy
 * - @ref DecayStrategy
 * - @ref PassthroughStrategy
 */


/**
 * @class ModelStrategy
 * @ingroup ParametersModule
 * @brief Abstract base class for different physics model strategies.
 */

class ModelStrategy {
public:
    /**
     * @brief Initializes model-specific parameters.
     * @param params Reference to Parameters object.
     * @return Set of block names that were absent during initialization.
     */
    virtual std::unordered_set<std::string> initializeParameters(class Parameters& params) = 0;

    /**
     * @brief Executes additional initialization tasks after main parameter loading.
     * @param params Reference to Parameters object.
     */
    virtual void postInitialization(Parameters& params) = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~ModelStrategy() = default;

    /**
     * @brief Sets the list of absent blocks for the model strategy.
     * @param _ Set of absent block names.
     */
    void add_absent_block(std::unordered_set<std::string> _) {absent_blocks = _;};

    /**
     * @brief Clears the list of absent blocks.
     */
    void remove_absent_block() {absent_blocks = std::unordered_set<std::string>();}
protected:
    std::unordered_set<std::string> absent_blocks;  ///< List of blocks missing from initialization.
};

/** 
 * @class SMModelStrategy
 * @ingroup ParametersModule
 * @brief Strategy for Standard Model parameters.
 */
class SMModelStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override;
};

/** 
 * @class BSMModelStrategy
 * @ingroup ParametersModule
 * @brief Strategy for Beyond Standard Model parameters.
 */
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

/** 
 * @class FlavorStrategy
 * @ingroup ParametersModule
 * @brief Strategy for flavor parameters (mesons mass, lifetime, etc.).
 */
class FlavorStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** 
 * @class GeneralModelStrategy
 * @ingroup ParametersModule
 * @brief Strategy for TODO ??.
 */
class GeneralModelStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** 
 * @class WilsonInputStrategy
 * @ingroup ParametersModule
 * @brief Strategy for Wilson parameters (and wilson itself).
 */
class WilsonInputStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** 
 * @class DecayStrategy
 * @ingroup ParametersModule
 * @brief Strategy for decay parameters.
 */
class DecayStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** 
 * @class ObservableStrategy
 * @ingroup ParametersModule
 * @brief Strategy for observables parameters (and observables themself).
 */
class ObservableStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/** 
 * @class PassthroughStrategy
 * @ingroup ParametersModule
 * @brief Strategy for passthrough parameters (not needed at runtime but must be in the output).
 */
class PassthroughStrategy : public ModelStrategy {
public:
std::unordered_set<std::string> initializeParameters(class Parameters& params) override;
    void postInitialization(Parameters& params) override {}
};

/**
 * @class Parameters
 * @ingroup ParametersModule
 * @brief Singleton class to manage parameter values, blocks, and model-specific strategies.
 */
class Parameters {
public:
    /**
     * @brief Retrieves the singleton instance for a given model.
     * @param id ParameterType (default: SM).
     * @return Shared pointer to the Parameters instance.
     */
    static std::shared_ptr<Parameters> GetInstance(ParameterType id = ParameterType::SM);

    /**
     * @brief Cleans up an instance of Parameters for a given model.
     * @param id ParameterType to remove.
     */
    void CleanupInstance(ParameterType id = ParameterType::SM);

    /**
     * @brief Checks if a parameter exists in a specified block.
     * @param block Block name.
     * @param pdgCode PDG code identifier.
     * @return True if parameter exists, false otherwise.
     */
    bool exist(const std::string& block, LhaID pdgCode);
    
    /**
     * @brief Retrieves a parameter value (operator syntax).
     * @param block Block name.
     * @param pdgCode PDG code identifier.
     * @return The corresponding parameter value.
     */
    scalar_t operator()(const std::string& block, LhaID pdgCode) const;

    /**
     * @brief Retrieves a shared pointer to a Parameter.
     * @param block Block name.
     * @param pdgCode PDG code identifier.
     * @return Shared pointer to the parameter.
     */
    std::shared_ptr<Parameter> get_parameter(const std::string& block, LhaID pdgCode);

    /**
     * @brief Sets a parameter value manually.
     * @param name Block name.
     * @param pdgCode PDG code.
     * @param value Value to assign.
     * @param force If true, forces the overwrite (default false).
     */
    void setBlockValue(const std::string& name, LhaID pdgCode, scalar_t value);

    /**
     * @brief Retrieves all parameter values for a given block.
     * @param blockName Block name.
     * @return Map of LHA IDs to parameter values.
     */
    std::map<LhaID, double> get_block_infos(std::string blockName);

    /**
     * @brief Retrieves the list of available parameter blocks.
     * @return Set of block names.
     */
    std::unordered_set<std::string> get_blocks_list();

    /**
     * @brief Changes the operational mode of a parameter (fixed/shiftable).
     * @param param_id Identifier of the parameter.
     * @param new_mode New mode to apply.
     */
    void changeParameterMode(const ParamId& param_id, ParameterMode new_mode);

    /**
     * @brief Applies a shift to a parameter value.
     * @param param_id Identifier of the parameter.
     * @param shift_value Value to shift.
     */
    void shiftParameter(const ParamId& param_id, scalar_t shift_value);

    /**
     * @brief Initializes parameter blocks for a given model type.
     * @param type Model type.
     * @return Set of missing blocks.
     */
    std::unordered_set<std::string> init_blocks(ParameterType type);

    /**
     * @brief Freezes an entire block (preventing parameter updates).
     * @param blockName Name of the block to freeze.
     */
    void freeze_block(const std::string& blockName);

    /**
     * @brief Unfreezes an entire block.
     * @param blockName Name of the block to unfreeze.
     */
    void unfreeze_block(const std::string& blockName);

    /**
     * @brief Freezes a specific parameter within a block.
     * @param blockName Name of the block.
     * @param id LHA ID of the parameter.
     */
    void freeze_param(const std::string& blockName, const LhaID& id);

    /**
     * @brief Unfreezes a specific parameter within a block.
     * @param blockName Name of the block.
     * @param id LHA ID of the parameter.
     */
    void unfreeze_param(const std::string& blockName, const LhaID& id);

    /**
     * @brief Destructor. Logs when a Parameters instance is destroyed.
     */
    ~Parameters() { LOG_DEBUG("Parameters at ", this); }

private:
    void claim_parameters(ParameterType type);
    
    /** 
     * @brief Private constructor for Parameters.
     * @param modelStrategy Strategy associated with this instance.
     */
    explicit Parameters(std::shared_ptr<ModelStrategy> modelStrategy);

    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;  ///< Static map of instances.
    std::shared_ptr<ModelStrategy> strategy;                                ///< Strategy object for this instance.
    std::shared_ptr<BlockAccessor> blockAccessor;                           ///< Accessor for blocks and parameters.

    /** @brief Factory friend. */
    friend class ParametersFactory;
    friend class DependentBlockManager;
};

/**
 * @class ParametersFactory
 * @ingroup ParametersModule
 * @brief Factory class for creating and managing Parameters instances.
 */
class ParametersFactory {
public:
    /**
     * @brief Retrieves the Parameters instance for a given model.
     * @param id ParameterType.
     * @return Shared pointer to Parameters.
     */
    static std::shared_ptr<Parameters> GetParameters(ParameterType id);

    /**
     * @brief Removes a Parameters instance for a given model.
     * @param id ParameterType.
     */
    static void removeParameters(ParameterType id);
private:
    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;  ///< Static map of instances.

    /**
     * @brief Creates the correct ModelStrategy for a given ParameterType.
     * @param id Model type.
     * @return Shared pointer to ModelStrategy.
     */
    static std::shared_ptr<ModelStrategy> createStrategy(ParameterType id);
};

#endif