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
#include "lha_reader.h"
#include "Interface.h"
#include "json_parser.h"
#include "QCDHelper.h"
#include "config.hpp"
#include <memory>
#include <ranges>
#include <algorithm>
#include <unordered_set>

struct ParameterBlockRepartition {
    static inline const std::map<ParameterType, std::vector<std::string>> BLOCKS {
        {ParameterType::SM, {"SMINPUTS", "MASS", "VCKMIN", "UPMNSIN", "GAUGE"}},
        {ParameterType::SUSY, {"MASS", "HMIX", "ALPHA", "MSOFT", "NMIX", "UMIX", "VMIX", "A0MIX", "H0MIX", "STOPMIX", "SBOTMIX", "STAUMIX", "AU", "AD", "AE", "YU", "YD", "YE"}},
        {ParameterType::THDM, {"MASS", "ALPHA", "LUCOUP", "DCOUPL", "LCOUPL"}},
        {ParameterType::FLAVOR, {"FMASS", "FLIFE", "FCONST", "FCONSTRATIO", "FBAG", "FPARAM"}},
        {ParameterType::WILSON, {"FWCOEF", "IMFWCOEF"}},
        {ParameterType::DECAY, {"B_Ks", "B_ll", "B_Xs", "B_Dlnu", "B_Dslnu"}},
        {ParameterType::OBSERVABLE, {"FOBS", "FOBSERR", "FOBSSM", "FDIPOLE"}},
        {ParameterType::PASSTHROUGH, {"MODSEL", "SPINFO", "FMODSEL", "FCINFO", "MINPAR", "EXTPAR"}}
    };

    static std::vector<std::string> filter_custom_blocks(const std::vector<std::string>& source);
};

struct ParametersAccessRights {
    static inline const std::map<std::string, std::unordered_set<int>> SM_RIGHTS {
        {"MASS", {1, 2, 3, 4, 11, 12, 13, 14, 15, 16, 21, 22, 24, 25}}, 
        {"GAUGE", {1, 2, 3}},
    };

    static inline const std::map<std::string, std::unordered_set<int>> THDM_RIGHTS {
        {"MASS", {25, 35, 36, 37}}, 
        {"GAUGE", {}},
    };

    static inline const std::map<std::string, std::unordered_set<int>> SUSY_RIGHTS {
        {"MASS", {25, 35, 36, 37, 
                  1000001, 1000002, 1000003, 1000004, 1000005, 1000006, 1000011, 1000012, 1000013, 1000014, 1000015, 1000016, 
                  2000001, 2000002, 2000003, 2000004, 2000005, 2000006, 2000011, 2000013, 2000015, 
                  1000021, 1000022, 1000023, 1000024, 1000025, 1000035, 1000037, 1000039}}, 
        {"GAUGE", {}},
    };
};

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
    virtual void initializeParameters(class Parameters& params) = 0;
    virtual ~ModelStrategy() = default;
};

// Concrete Strategy Classes
/** @class SMModelStrategy @brief Strategy for the Standard Model. */
class SMModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

/** @class SUSYModelStrategy @brief Strategy for SUSY models. */
class SUSYModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

/** @class THDMModelStrategy @brief Strategy for Two-Higgs-Doublet Models. */
class THDMModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

/** @class FlavorStrategy @brief Strategy for flavor physics parameters. */
class FlavorStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

/** @class GeneralModelStrategy @brief General model strategy for parameter initialization. */
class GeneralModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

/** @class WilsonInputStrategy @brief Strategy for Wilson coefficient inputs. */
class WilsonInputStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

/** @class FormFactorStrategy @brief Strategy for form factor parameters. */
class FormFactorStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
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
     * @brief Determines the parameter type based on a block and PDG code.
     * @param block The name of the parameter block.
     * @param pdgCode The PDG code associated with the parameter.
     * @return The corresponding ParameterType.
     */
    static ParameterType GetType(const std::string& block, LhaID pdgCode);

    /**
     * @brief Retrieves a parameter value given its type, block, and code.
     * @param type The parameter type.
     * @param block The name of the block containing the parameter.
     * @param code The parameter code.
     * @return The retrieved parameter value.
     */
    static double Get(ParameterType type, const std::string& block, LhaID code);

    /**
     * @brief Retrieves a parameter value given its unique identifier.
     * @param id The unique identifier of the parameter.
     * @return The parameter value.
     */
    static double Get(ParamId id);

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
    double operator()(const std::string& block, LhaID pdgCode);

    /**
     * @brief Adds a new block of parameters to the collection.
     * @param name The name of the new block.
     * @param block A shared pointer to the block.
     */
    void addBlock(const std::string& name, std::shared_ptr<Block> block);

    /**
     * @brief Adds a new block dependant from another of parameters to the collection.
     * @param name The name of the new block.
     * @param block A shared pointer to the block.
     * @param source_block The name of the block this one depend one.
     * @param recalculateFunc The function that describe the dependancy.
     */
    void addDependantBlock(const std::string& name, std::shared_ptr<DependentBlock>& block,const std::string& source_block, std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateFunc);

    /**
     * @brief Sets a parameter value within a specified block.
     * @param name The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to assign.
     * @param force If true, forces the update.
     */
    void setBlockValue(const std::string& name, LhaID pdgCode, double value, bool force = false);

    /**
     * @brief Retrieves all parameter values from a specified block.
     * @param blockName The name of the block.
     * @return A map of PDG codes to parameter values.
     */
    std::map<LhaID, double> get_block_infos(std::string blockName);

    /**
     * @brief Retrieves a list of all available parameter blocks.
     * @return A vector of block names.
     */
    std::vector<std::string> get_blocks_list();

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

    /**
     * @brief Retrieves a specific entry from the CKM matrix.
     * @param idx The index of the matrix entry.
     * @return The corresponding complex CKM entry.
     */
    static complex_t get_c_CKM_entry(LhaID idx);

    void init_blocks(ParameterType type, std::unordered_set<std::string> excluded_dependant = {});

    /**
     * @brief Destructor that logs the destruction of a Parameters instance.
     */
    ~Parameters() { LOG_DEBUG("Parameters at ", this); }

private:
    /** @brief Private constructor for singleton pattern. */
    explicit Parameters(std::shared_ptr<ModelStrategy> modelStrategy);

    /** @brief Map of model types to their corresponding Parameters instances. */
    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;

    /** @brief Strategy used for parameter management. */
    std::shared_ptr<ModelStrategy> strategy;

    std::shared_ptr<BlockAccessor> blockAccessor;

    /** @brief Factory friend. */
    friend class ParametersFactory;
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