#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <memory>
#include <ranges>
#include <algorithm>
#include <unordered_set>
#include <functional>

#include "ParameterRouter.h"
#include "BlockAccessor.h"
#include "Interface.h"
#include "QCDHelper.h"
#include "EWHelper.h"
#include "config.hpp"
#include "SourcesView.h"
#include "DependentParameter.h"

// Forward declaration for strategy interface
class Parameters;

/**
 * @file Parameters.h
 * @brief Model-dependent parameter repository and initialization strategies.
 *
 * This header defines the main runtime parameter layer used throughout HyperISO.
 *
 * It contains:
 * - @ref ModelStrategy, the abstract interface used to initialize one parameter domain,
 * - concrete strategies for each @ref ParameterType
 *   (SM, BSM, FLAVOR, WILSON, DECAY, OBSERVABLE, PASSTHROUGH),
 * - @ref Parameters, the main block-based repository for one parameter domain,
 * - @ref ParametersFactory, the singleton-style factory ensuring one repository
 *   per @ref ParameterType.
 *
 * ## General architecture
 *
 * A @ref Parameters instance owns a @ref BlockAccessor containing all blocks
 * for one parameter namespace (SM, BSM, WILSON, ...). The actual content is
 * initialized through a model-specific @ref ModelStrategy:
 *
 * @code
 *   Parameters::GetInstance(ParameterType::SM)
 *       -> ParametersFactory::GetParameters(ParameterType::SM)
 *       -> createStrategy(ParameterType::SM)
 *       -> strategy->initializeParameters(...)
 *       -> strategy->postInitialization(...)
 * @endcode
 *
 * Blocks are usually extracted from the @ref MemoryManager input cache, then
 * enriched with dependent blocks / derived quantities during post-initialization.
 *
 * ## Responsibilities
 *
 * @ref Parameters provides:
 * - block and parameter lookup,
 * - manual value updates,
 * - parameter/block freezing,
 * - dependency detachment / reattachment,
 * - access to the underlying @ref BlockAccessor for advanced operations.
 *
 * @see Parameters
 * @see ParametersFactory
 * @see ModelStrategy
 * @see BlockAccessor
 * @see MemoryManager
 */

/**
 * @defgroup ParametersModule Parameters and model strategies
 * @brief Management of physics parameters and model-specific initialization logic.
 *
 * This module groups the parameter repositories and the initialization strategies
 * used to populate them from the input cache and derived helper computations.
 */

/**
 * @class ModelStrategy
 * @ingroup ParametersModule
 * @brief Abstract initialization strategy for one parameter namespace.
 *
 * A ModelStrategy is responsible for two phases:
 * - @ref initializeParameters: load the initial set of blocks for a given
 *   parameter type,
 * - @ref postInitialization: build extra derived blocks / helper quantities
 *   after the raw blocks have been loaded.
 *
 * It also stores the set of blocks that were absent during initialization.
 */
class ModelStrategy {
public:
    /**
     * @brief Initializes the block set for the target parameter repository.
     *
     * Typical implementations call @ref Parameters::init_blocks with the
     * appropriate @ref ParameterType and return the set of missing blocks.
     *
     * @param params Parameter repository being initialized.
     * @return Set of block names that were expected but absent from the input.
     */
    virtual std::unordered_set<BlockName> initializeParameters(class Parameters& params) = 0;

    /**
     * @brief Performs additional initialization after raw blocks are loaded.
     *
     * This hook is typically used to:
     * - initialize helper modules,
     * - create derived / dependent blocks,
     * - add fallback blocks when some inputs are absent,
     * - synchronize scales or helper quantities.
     *
     * @param params Parameter repository being post-initialized.
     */
    virtual void postInitialization(Parameters& params) = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~ModelStrategy() = default;

    /**
     * @brief Stores the set of blocks absent during initialization.
     * @param _ Set of missing block names.
     */
    void add_absent_block(std::unordered_set<BlockName> _) {absent_blocks = _;};

    /**
     * @brief Clears the list of absent blocks.
     */
    void remove_absent_block() {absent_blocks = std::unordered_set<BlockName>();}
protected:
    std::unordered_set<BlockName> absent_blocks;  /// Set of blocks that were expected but not present in the input cache.
};

/**
 * @class SMModelStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for Standard Model parameters.
 *
 * Raw SM blocks are loaded from the input cache through
 * @ref Parameters::init_blocks(ParameterType::SM).
 *
 * Post-initialization currently:
 * - initializes @ref QCDHelper and @ref EWHelper,
 * - creates a derived `VCKM` block from `VCKMIN` if `VCKM` is absent,
 * - ensures the Wilson parameter repository exists,
 * - creates a derived `MASS_EW_SCALE` block using `EW_SCALE` and QCD running masses.
 */
class SMModelStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     */
    void postInitialization(Parameters& params) override;
};

/**
 * @class BSMModelStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for BSM parameters.
 *
 * Raw BSM blocks are loaded from the input cache through
 * @ref Parameters::init_blocks(ParameterType::BSM).
 *
 * In the current implementation, post-initialization contains SUSY-specific
 * fallback logic: if some squark-mixing blocks are absent, zero-filled
 * dependent 7x7 blocks are created.
 */
class BSMModelStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     */
    void postInitialization(Parameters& params) override;
};

/**
 * @class FlavorStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for flavor-sector input blocks.
 *
 * This strategy loads the FLAVOR parameter namespace from the input cache.
 * It currently has no extra post-initialization step.
 */
class FlavorStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     *
     * Current implementation does nothing.
     */
    void postInitialization(Parameters&) override {}
};

/**
 * @class GeneralModelStrategy
 * @ingroup ParametersModule
 * @brief Placeholder / generic strategy class.
 *
 * This strategy is declared in the interface but its behavior is not shown in
 * the implementation snippet provided here. It appears intended as a generic
 * or catch-all strategy for model-specific blocks outside the standard sets.
 *
 * @warning If kept in the public API, its implementation should be documented
 *          and provided consistently.
 */
class GeneralModelStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     *
     * Current implementation does nothing.
     */
    void postInitialization(Parameters&) override {}
};

/**
 * @class WilsonInputStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for Wilson-coefficient input blocks.
 *
 * This strategy loads the WILSON namespace from the input cache.
 *
 * In the current implementation, post-initialization additionally synchronizes
 * the helper block `EW_SCALE` with the scale carried by the `FWCOEF` block,
 * but only if the configuration flag
 * @ref ExternalFlag::HAS_WILSON_INPUT is enabled.
 */
class WilsonInputStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     */
    void postInitialization(Parameters& params) override;
};

/**
 * @class DecayStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for decay-related parameter blocks.
 *
 * This strategy loads the DECAY namespace from the input cache and performs
 * no extra post-initialization step in the current implementation.
 */
class DecayStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     *
     * Current implementation does nothing.
     */
    void postInitialization(Parameters&) override {}
};

/**
 * @class ObservableStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for observable-related blocks.
 *
 * This strategy loads the OBSERVABLE namespace from the input cache and
 * currently performs no additional post-initialization.
 */
class ObservableStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     *
     * Current implementation does nothing.
     */
    void postInitialization(Parameters&) override {}
};

/**
 * @class PassthroughStrategy
 * @ingroup ParametersModule
 * @brief Strategy responsible for passthrough / output-preserved blocks.
 *
 * This namespace is intended for blocks that are not actively used by the
 * runtime logic but still need to be carried through the framework and/or
 * reproduced in outputs.
 */
class PassthroughStrategy : public ModelStrategy {
public:
    /**
     * @copydoc ModelStrategy::initializeParameters
     */
    std::unordered_set<BlockName> initializeParameters(class Parameters& params) override;

    /**
     * @copydoc ModelStrategy::postInitialization
     *
     * Current implementation does nothing.
     */
    void postInitialization(Parameters&) override {}
};

/**
 * @class Parameters
 * @ingroup ParametersModule
 * @brief Block-based parameter repository for one @ref ParameterType namespace.
 *
 * A Parameters instance owns the complete block set for one parameter domain
 * (SM, BSM, FLAVOR, WILSON, ...), together with the strategy that initialized it.
 *
 * Responsibilities include:
 * - accessing scalar values and full Parameter objects,
 * - mutating values,
 * - freezing/unfreezing blocks and parameters,
 * - detaching/reattaching dependencies,
 * - listing blocks,
 * - exposing the underlying @ref BlockAccessor for advanced workflows.
 *
 * Instances are obtained through @ref GetInstance and are managed through
 * @ref ParametersFactory, effectively giving one repository per @ref ParameterType.
 */
class Parameters {
public:
    /**
     * @brief Returns the singleton-like repository for a given parameter type.
     *
     * This first checks whether the requested parameter type is enabled in the
     * current @ref MemoryManager cache, then delegates creation/retrieval to
     * @ref ParametersFactory.
     *
     * @param id Parameter namespace to retrieve (default: SM).
     * @return Shared pointer to the repository.
     *
     * @throws Via logging/error handling if the parameter type is not allowed
     *         by the current memory cache.
     */
    static std::shared_ptr<Parameters> GetInstance(ParameterType id = ParameterType::SM);

    /**
     * @brief Removes one repository instance from the factory cache.
     *
     * This is a thin wrapper around @ref ParametersFactory::removeParameters.
     *
     * @param id Parameter namespace to remove.
     */
    void CleanupInstance(ParameterType id = ParameterType::SM);

    /**
     * @brief Checks whether a parameter exists in a given block.
     *
     * Current implementation:
     * - returns false if the block does not exist,
     * - otherwise delegates to `Block::contains(id)`.
     *
     * @param block Block name.
     * @param pdgCode Parameter identifier inside the block.
     * @return True if the parameter exists, false otherwise.
     */
    bool exist(const BlockName& block, LhaID pdgCode);
    
    /**
     * @brief Reads the current scalar value of a parameter.
     *
     * This is a convenience wrapper over the underlying @ref BlockAccessor.
     *
     * @param block Block name.
     * @param pdgCode Parameter identifier.
     * @return Current scalar value.
     */
    scalar_t operator()(const BlockName& block, LhaID pdgCode) const;

    /**
     * @brief Returns the full @ref Parameter object for one entry.
     *
     * @param block Block name.
     * @param pdgCode Parameter identifier.
     * @return Shared pointer to the stored parameter.
     */
    std::shared_ptr<Parameter> get_parameter(const BlockName& block, LhaID pdgCode);

    /**
     * @brief Sets or creates a scalar value inside a block.
     *
     * Delegates to `BlockAccessor::setValue(...)`.
     *
     * @param name Block name.
     * @param pdgCode Parameter identifier.
     * @param value Value to set.
     */
    void setBlockValue(const BlockName& name, LhaID pdgCode, scalar_t value);

    /**
     * @brief Returns all scalar values stored in one block.
     *
     * @param blockName Block name.
     * @return Map of `LhaID -> scalar_t`.
     */
    std::map<LhaID, scalar_t> get_block_infos(BlockName blockName);

    /**
     * @brief Returns the scale attached to a block.
     *
     * This directly delegates to `Block::get_scale()`.
     *
     * @param blockName Block name.
     * @return Block scale.
     *
     * @throws If the block does not exist or the block has no scale set.
     */
    double get_block_scale(BlockName blockName) const;
    

    /**
     * @brief Returns the list of blocks currently held in this repository.
     *
     * @return Set of block names.
     */
    std::unordered_set<BlockName> get_blocks_list();

    /**
     * @brief Checks whether one block is a dependent block.
     *
     * This delegates to @ref BlockAccessor::is_dependent_block on the repository
     * owned by this @ref Parameters instance.
     *
     * @param blockName Block name or alias to inspect.
     * @return True if the block is a @ref DependentBlock, false otherwise.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    bool is_dependent_block(const BlockName& blockName) const;

    /**
     * @brief Returns the direct source blocks of one block.
     *
     * Direct source blocks are the upstream blocks immediately used by the
     * inspected block. Plain blocks return an empty list.
     *
     * @param blockName Block name or alias to inspect.
     * @return Sorted list of direct upstream/source block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_source_blocks(const BlockName& blockName) const;

    /**
     * @brief Returns the direct blocks depending on one block.
     *
     * This follows direct block observers and is the downstream counterpart of
     * @ref get_source_blocks.
     *
     * @param blockName Block name or alias to inspect.
     * @return Sorted list of direct downstream/dependent block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_dependent_blocks(const BlockName& blockName) const;

    /**
     * @brief Returns all transitive source blocks of one block.
     *
     * This recursively walks upstream through dependent-block source links.
     *
     * @param blockName Block name or alias to inspect.
     * @return Sorted list of all upstream/source block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_all_source_blocks(const BlockName& blockName) const;

    /**
     * @brief Returns all transitive blocks depending on one block.
     *
     * This recursively walks downstream through block observers.
     *
     * @param blockName Block name or alias to inspect.
     * @return Sorted list of all downstream/dependent block names.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::vector<std::string> get_all_dependent_blocks(const BlockName& blockName) const;
    
    /**
     * @brief Changes the mode of one parameter.
     *
     * @warning This is currently not implemented in the provided `.cpp`
     *          (the method body is effectively a no-op).
     *
     * @param param_id Parameter identifier.
     * @param new_mode New parameter mode.
     */
    void changeParameterMode(const ParamId& param_id, ParameterMode new_mode);

    /**
     * @brief Applies an additive shift to a parameter by rewriting its value.
     *
     * Current implementation simply does:
     * @code
     * value <- current_value + shift_value
     * @endcode
     * through the accessor, rather than using `Parameter::set_shift()`.
     *
     * @param param_id Parameter identifier.
     * @param shift_value Additive shift.
     */
    void shiftParameter(const ParamId& param_id, scalar_t shift_value);

    /**
     * @brief Loads the raw blocks for one parameter namespace from MemoryManager.
     *
     * This method:
     * - compares the expected blocks from `ParameterBlockRepartition::BLOCKS`
     *   against the input cache,
     * - separates present and missing blocks,
     * - applies a few policy adjustments (e.g. FWCOEF handling),
     * - extracts the present blocks through `MemoryManager::extract_blocks(...)`,
     * - tags all extracted parameters with the correct owner type.
     *
     * @param type Parameter namespace being initialized.
     * @return Set of missing blocks.
     */
    std::unordered_set<BlockName> init_blocks(ParameterType type);

    /**
     * @brief Freezes one whole block.
     *
     * Delegates to `Block::freeze()` after checking existence.
     *
     * @param blockName Block to freeze.
     */
    void freeze_block(const BlockName& blockName);

    /**
     * @brief Freezes one whole block.
     *
     * Delegates to `Block::freeze()` after checking existence.
     *
     * @param blockName Block to freeze.
     */
    void unfreeze_block(const BlockName& blockName);

    /**
     * @brief Freezes one parameter inside a block.
     *
     * @param blockName Block containing the parameter.
     * @param id Parameter identifier.
     */
    void freeze_param(const BlockName& blockName, const LhaID& id);

    /**
     * @brief Unfreezes one parameter inside a block.
     *
     * @param blockName Block containing the parameter.
     * @param id Parameter identifier.
     */
    void unfreeze_param(const BlockName& blockName, const LhaID& id);

    /**
     * @brief Detaches one dependent block from its upstream dependency graph.
     *
     * Delegates to `BlockAccessor::detach_block(...)`.
     *
     * @param blockName Block to detach.
     */
    void detach_block(const BlockName& blockName);

    /**
     * @brief Reattaches one previously detached dependent block.
     *
     * Delegates to `BlockAccessor::reattach_block(...)`.
     *
     * @param blockName Block to reattach.
     */
    void reattach_block(const BlockName& blockName);

    /**
     * @brief Detaches one dependent parameter from its upstream dependencies.
     *
     * Delegates to `BlockAccessor::detach_parameter(...)`.
     *
     * @param blockName Block containing the parameter.
     * @param id Parameter identifier.
     */
    void detach_param(const BlockName& blockName, const LhaID& id);

    /**
     * @brief Reattaches one previously detached dependent parameter.
     *
     * Delegates to `BlockAccessor::reattach_parameter(...)`.
     *
     * @param blockName Block containing the parameter.
     * @param id Parameter identifier.
     */
    void reattach_param(const BlockName& blockName, const LhaID& id);

    /**
     * @brief Returns the underlying block accessor.
     *
     * This is mainly intended for advanced / optimization-oriented workflows.
     *
     * @return Shared pointer to the block accessor.
     */
    std::shared_ptr<BlockAccessor> get_block_accessor() { return this->blockAccessor;}

    /**
     * @brief Prints one block to stdout.
     *
     * Internally streams the block object resolved by name through
     * the block accessor.
     *
     * @param blockname Name of the block to print.
     */
    void print_block(const std::string blockname);
    
    /**
     * @brief Stream output operator for the whole repository.
     *
     * Prints the underlying @ref BlockAccessor.
     *
     * @param os Output stream.
     * @param instance Repository instance to print.
     * @return Output stream.
     */
    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<Parameters>);
    
    /**
     * @brief Destructor.
     *
     * Only logs destruction in the current implementation.
     */
    ~Parameters() { LOG_DEBUG("Parameters at ", this); }

private:
    /**
     * @brief Tags all currently stored parameters with the given owner type.
     *
     * This iterates through all blocks and calls `Block::set_owner(type)`.
     *
     * @param type Parameter owner type to assign.
     */
    void claim_parameters(ParameterType type);
    
    /**
     * @brief Constructs a repository using a model strategy.
     *
     * The constructor immediately calls `strategy->initializeParameters(*this)`
     * and stores the returned missing-block set into the strategy.
     *
     * @param modelStrategy Strategy used to initialize this repository.
     */
    explicit Parameters(std::shared_ptr<ModelStrategy> modelStrategy);

    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;  /// Static cache of repository instances (declared here, factory-managed in practice).
    std::shared_ptr<ModelStrategy> strategy;                                /// Initialization strategy associated with this repository.
    std::shared_ptr<BlockAccessor> blockAccessor;                           /// Underlying storage/access layer for all parameter blocks.

    /** @brief Factory friend. */
    friend class ParametersFactory;                                          /// Factory is allowed to construct/remove repositories.
    friend class DependentBlockManager;                                      /// Dependent block manager may access internals when wiring derived blocks.
};

/**
 * @class ParametersFactory
 * @ingroup ParametersModule
 * @brief Factory/cache for @ref Parameters instances.
 *
 * This class ensures that at most one repository exists per @ref ParameterType.
 * It is responsible for:
 * - creating the matching @ref ModelStrategy,
 * - constructing the repository,
 * - calling strategy post-initialization,
 * - caching the result for later reuse.
 */
class ParametersFactory {
public:
    /**
     * @brief Returns the cached repository for one parameter type.
     *
     * If it does not exist yet:
     * - creates the appropriate strategy,
     * - constructs the repository,
     * - runs post-initialization,
     * - stores it in the static cache.
     *
     * @param id Parameter namespace to retrieve.
     * @return Shared pointer to the corresponding repository.
     */
    static std::shared_ptr<Parameters> GetParameters(ParameterType id);

    /**
     * @brief Creates a fresh uncached repository for a runtime context.
     *
     * The returned repository is initialized and post-initialized exactly like
     * the global singleton repository, but it is not stored in the process-wide
     * factory cache.
     */
    static std::shared_ptr<Parameters> CreateUncached(ParameterType id);

    /**
     * @brief Creates a fresh uncached repository and registers it before post-initialization.
     *
     * Some post-initialization steps create dependent blocks through
     * DependentBlockManager, which can recursively call Parameters::GetInstance
     * for the repository currently being built. Runtime-local contexts need the
     * partially constructed repository to be visible before those dependent
     * blocks are attached. The callback is invoked after the base repository is
     * initialized but before ModelStrategy::postInitialization runs.
     */
    static std::shared_ptr<Parameters> CreateUncachedRegistered(
        ParameterType id,
        const std::function<void(const std::shared_ptr<Parameters>&)>& register_before_postinit
    );

    /**
     * @brief Removes one cached repository.
     *
     * @param id Parameter namespace to remove.
     *
     * @throws Via logging/error handling if the repository does not exist.
     */
    static void removeParameters(ParameterType id);

    static void clear();
private:
    static std::map<ParameterType, std::shared_ptr<Parameters>> instances;  /// Static cache of repositories indexed by parameter type.

    /**
     * @brief Creates the concrete strategy associated with a parameter type.
     *
     * Current mapping:
     * - SM          -> @ref SMModelStrategy
     * - BSM         -> @ref BSMModelStrategy
     * - FLAVOR      -> @ref FlavorStrategy
     * - WILSON      -> @ref WilsonInputStrategy
     * - DECAY       -> @ref DecayStrategy
     * - OBSERVABLE  -> @ref ObservableStrategy
     * - PASSTHROUGH -> @ref PassthroughStrategy
     *
     * @param id Parameter namespace.
     * @return Shared pointer to the matching strategy.
     *
     * @throws std::invalid_argument for unsupported parameter types.
     */
    static std::shared_ptr<ModelStrategy> createStrategy(ParameterType id);
};

#endif