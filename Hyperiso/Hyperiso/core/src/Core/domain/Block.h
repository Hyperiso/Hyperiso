#ifndef BLOCK_H
#define BLOCK_H

#include <array>
#include <functional>

#include "Include.h"
#include "Parameter.h"
#include "IStorage.h"

/**
 * @file Block.h
 * @brief Defines classes for managing parameter blocks and their dependencies.
 *
 * This file declares a hierarchy of classes used to store, manage, and update
 * parameter values across different parameter instances. It includes:
 *
 * - Block:      a container of parameters identified by LHA IDs.
 * - DependentBlock: a Block whose content is computed from one or more source Blocks.
 *
 * Blocks support ownership tagging, global scales, a simple observer mechanism
 * (between blocks and between parameters), and dependency propagation with
 * freeze/unfreeze semantics.
 */

class BlockSrc;
class Block;
class DependentBlock;

/**
 * @typedef DepUpdateFunc
 * @brief Function type used to update a DependentBlock from its sources.
 *
 * The function receives:
 *  - a BlockSrc view over the source blocks,
 *  - a shared pointer to the DependentBlock being updated.
 *
 * It is responsible for recomputing the content of the dependent block
 * (e.g. filling parameters, updating values, etc.).
 */
typedef std::function<void(const BlockSrc&, std::shared_ptr<DependentBlock>)> DepUpdateFunc;


/**
 * @class Block
 * @brief Container for a collection of parameters identified by LHA IDs.
 *
 * A Block:
 *  - stores parameters in a map keyed by @ref LhaID,
 *  - implements the @ref IStorage interface,
 *  - can notify other blocks (observers) when its content changes,
 *  - can propagate freeze/unfreeze/update calls down to its parameters,
 *  - can carry an optional global renormalization scale.
 *
 * Blocks are typically grouped by @ref BlockName (e.g. an SLHA or FLHA block).
 */
class Block : public IStorage<LhaID, Parameter>, public std::enable_shared_from_this<Block> {
public:
    /// Name of the block (e.g. "SMINPUTS", "MASS", "FWCOEF"...)
    BlockName blockname {""};

    /// Default constructor (empty block with empty name).
    Block() = default;

    /**
     * @brief Copy-like constructor from another block.
     *
     * Copies:
     *  - items map,
     *  - block name,
     *  - scale (if present),
     * and re-assigns the owner_block pointer of each parameter to this new block.
     *
     * @param other Shared pointer to the source block.
     */
    Block(std::shared_ptr<Block> other);

    /**
     * @brief Returns the block name.
     */
    BlockName get_name() const { return this->blockname; }

    /**
     * @brief Stores a new parameter in the block.
     *
     * If the key already exists, the parameter is not overwritten and a debug
     * message is logged.
     *
     * @param id    The LHA ID of the parameter.
     * @param param Shared pointer to the Parameter to store.
     */
    void store(const LhaID& id, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Overwrites the numerical payload of an existing entry.
     *
     * This does NOT replace the stored shared_ptr. Instead it calls
     * Parameter::overwrite_payload_from() on the existing parameter, i.e. it copies
     * expected/std/shift/mode/scale/binning, but keeps identity (pointer), id,
     * observers, and owner_block of the destination object.
     *
     * After overwriting:
     *  - notifies parameter observers via dst->notifyObservers(),
     *  - then notifies block observers via Block::notifyObservers().
     */
    void assign(const LhaID& key, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Assigns a new numerical value to an existing parameter.
     *
     * Implementation note (current behavior):
     *  - calls Parameter::set_expected(value) which already triggers parameter-level
     *    observer notifications (and may notify the owning block via owner_block),
     *  - then explicitly calls Parameter::notifyObservers() again,
     *  - then notifies observer blocks via Block::notifyObservers().
     *
     * @param key   The LHA ID of the parameter.
     * @param value The new central value to assign.
     */
    void assign(const LhaID& key, scalar_t value);

    /**
     * @brief Stores or assigns a parameter depending on whether it already exists.
     *
     * Note: when storing a new parameter, this block will set the parameter's
     * owner_block (if a valid self weak pointer is available).
     */
    void store_or_assign(const LhaID& key, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Checks if a parameter exists in the block.
     * @param key The LHA ID of the parameter.
     * @return True if the parameter exists, false otherwise.
     */
    bool contains(const LhaID& key) const override;

    /**
     * @brief Retrieves a parameter by ID (ensures the block is up-to-date first).
     *
     * Calls ensure_up_to_date() before lookup.
     * If the key is missing, an error is logged; depending on the logging policy,
     * the program may abort/throw. The subsequent items.at(id) would also throw
     * std::out_of_range if reached.
     */
    std::shared_ptr<Parameter> retrieve(const LhaID& id) override;

    /**
     * @brief Removes a parameter and clears its dependency tree.
     *
     * The removed parameter has its @ref Parameter::clear_below() called,
     * then it is erased from the block. Afterwards, all dependent blocks
     * (observers) are destroyed recursively.
     *
     * @param key The LHA ID of the parameter to remove.
     */
    void remove(const LhaID& key) override;

    /**
     * @brief Retrieves the set of all parameter IDs stored in the block.
     * @return Set of LhaID keys.
     */
    std::unordered_set<LhaID> getAllIDs();

    /**
     * @brief Provides access to the internal map of stored parameters.
     *
     * Used mainly for iteration, debugging, and printing.
     *
     * @return Const reference to the items map.
     */
    const std::map<LhaID, std::shared_ptr<Parameter>>& getItems();

    /**
     * @brief Sets the owner type (ParameterType) for all parameters in the block.
     *
     * Internally calls @ref Parameter::set_owner() on each entry.
     *
     * @param type The ParameterType to set as owner.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Adds an observer block that will be notified on updates.
     *
     * The observer's @ref Block::update() method is called by @ref notifyObservers().
     *
     * @param observer Shared pointer to the observer block.
     */
    void addObserver(std::shared_ptr<Block> observer);

    /**
     * @brief Removes an observer block if present.
     * @param observer Shared pointer to the observer block.
     */
    void removeObserver(std::shared_ptr<Block> observer);

    /**
     * @brief Notifies all observer blocks that this block has been updated.
     *
     * For each non-null observer, calls @ref Block::update().
     * Null observers are cleaned up from the list.
     */
    void notifyObservers();

    /**
     * @brief Retrieves the list of observer blocks.
     * @return Vector of shared pointers to observer blocks.
     */
    std::vector<std::shared_ptr<Block>> getObservers() const;

    /**
     * @brief Virtual hook to initialize the block.
     *
     * Default behavior: do nothing. Derived classes may override to
     * perform initial filling, dependency wiring, etc.
     */
    virtual void init() {};
    
    /**
     * @brief Updates the block content.
     *
     * Default behavior: calls @ref Parameter::update() on all stored parameters.
     * DependentBlock overrides this to recompute its content from source blocks.
     */
    virtual void update();

    /**
     * @brief Destroys this block and its dependency subtree.
     *
     * Current order:
     *  - takes ownership of the observer list (dependents),
     *  - calls clear_below() on all parameters (snapshot-based),
     *  - clears the local items map,
     *  - recursively calls destroy() on dependent blocks.
     *
     * @note This function prints to std::cout (currently).
     */
    virtual void destroy();
    
    /**
     * @brief Freezes all parameters in the block.
     *
     * Default behavior: calls @ref Parameter::freeze() on each stored parameter.
     * Derived classes may extend this behavior.
     */
    virtual void freeze();

    /**
     * @brief Unfreezes all parameters in the block.
     *
     * Default behavior: calls @ref Parameter::unfreeze() on each stored parameter.
     */
    virtual void unfreeze();

    /**
     * @brief Copies the content and metadata from another block.
     *
     * Copies:
     *  - items map (shared_ptr<Parameter>),
     *  - blockname,
     *  - scale (if present).
     *
     * Also updates each parameter's owner_block to this block.
     *
     * @param other Shared pointer to the source block.
     */
    void copy(std::shared_ptr<Block> other);

    /**
     * @brief Clears dependencies above all parameters in this block.
     *
     * By default, calls @ref Parameter::clear_above() on each stored parameter.
     */
    virtual void clear_above();

    /**
     * @brief Clears dependencies below all parameters in this block.
     *
     * Implementation detail:
     *  - builds a snapshot vector of the current parameter pointers,
     *    then calls clear_below() on each pointer.
     * This avoids iterator invalidation if dependency clearing mutates structures.
     */
    virtual void clear_below();

    /**
     * @brief Checks whether the block has a defined global scale.
     * @return True if the scale is set, false otherwise.
     */
    bool has_scale();

    /**
     * @brief Sets the global scale of the block.
     *
     * @param scale Scale value to set.
     * @throws If a scale is already set, logs a logic error.
     */
    void set_scale(double scale);

    /**
     * @brief Retrieves the global scale of the block.
     *
     * @return The scale value.
     * @throws If no scale is set, logs a logic error.
     */
    double get_scale();

    /**
     * @brief Erases a parameter without touching its dependencies.
     *
     * Unlike @ref remove(), this does not call clear_below() on the
     * parameter nor destroy dependent blocks. It strictly erases the
     * entry in the local items map.
     *
     * @param id LhaID of the parameter to erase.
     */
    void erase_local(const LhaID& id);

    /**
     * @brief Returns the source blocks of this block.
     *
     * For a plain Block, there are no source blocks and the default
     * implementation returns an empty map. DependentBlock overrides this.
     *
     * @return Map of source block names to their shared pointers.
     */
    virtual std::unordered_map<std::string, std::shared_ptr<Block>> get_source_blocks() const;

    //TODO 
    virtual void ensure_up_to_date() {}

    void bind_self(std::shared_ptr<Block> self) {
        self_ = self;

        std::weak_ptr<Block> w(self);
        for (auto& [_, p] : items) {
            if (p) p->set_owner_block(w);
        }
    }

    /**
     * @brief Destructor.
     *
     * Does not automatically destroy dependencies; explicit @ref destroy()
     * is used for that. Kept virtual through inheritance.
     */
    ~Block() {}

    /**
     * @brief Pretty-print utility for a Block.
     *
     * Outputs:
     * @code
     * Block NAME:
     *     id: value
     *     ...
     * @endcode
     */
    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<Block>);
protected:
    std::vector<std::shared_ptr<Block>> observers;      ///< List of observer blocks.
    std::map<LhaID, std::shared_ptr<Parameter>> items;  ///< Map of parameters stored in this block.
    std::optional<double> scale;                        ///< Optional global scale associated with this block.

    std::weak_ptr<Block> self_;

    std::weak_ptr<Block> self_weak() {
        if (auto s = self_.lock()) return self_;
        try { return shared_from_this(); } catch (...) { return {}; }
    }
};

/**
 * @class DependentBlock
 * @brief Block whose content depends on one or more source Blocks.
 *
 * A DependentBlock:
 *  - keeps a map of source blocks (@ref sourceBlocks),
 *  - registers itself as observer of its sources,
 *  - recomputes its own content via a user-supplied @ref DepUpdateFunc,
 *  - supports freeze/unfreeze logic with deferred updates.
 *
 * It is typically used to represent derived quantities (e.g. running parameters,
 * translated bases, combined uncertainties, etc.) built from primary blocks.
 */
class DependentBlock : public Block {
public:
    /**
     * @brief Constructs a DependentBlock.
     *
     * @param sources        Map from block names to shared pointers of source Blocks.
     * @param recalculateFunc Function used to recompute this block from its sources.
     */
    explicit DependentBlock(const std::unordered_map<std::string, std::shared_ptr<Block>>& sources, DepUpdateFunc recalculateFunc) 
        : sourceBlocks(std::move(sources)), recalculateLambda(std::move(recalculateFunc)), frozen(false) {}

    /**
     * @brief Checks if this block depends on a specific source block.
     * @param blockName Name of the source block.
     * @return True if a non-null source block with this name exists, false otherwise.
     */
    bool dependsOn(const std::string& blockName);

    /**
     * @brief Initializes dependency tracking.
     *
     * Must be called after construction. Registers this DependentBlock as an
     * observer of all source blocks.
     */
    void init() override;

    /**
     * @brief Marks the block as dirty (lazy update).
     *
     * This override does not recompute immediately.
     *  - If frozen: sets a pending flag and returns.
     *  - Otherwise: propagates dirtiness to dependent DependentBlocks via mark_dirty().
     *
     * Actual recomputation is performed on demand in ensure_up_to_date_impl(),
     * typically triggered by calling Block::retrieve/contains/getItems/getAllIDs.
     */
    void update() override;

    /**
     * @brief Freezes the block (no updates will be performed).
     *
     * While frozen, calls to @ref update() will only set a flag indicating
     * that an update is pending.
     */
    void freeze() override;

    /**
     * @brief Unfreezes the block and triggers an update if one was pending.
     */
    void unfreeze() override;

    /**
     * @brief Destructor.
     *
     * Unregisters this block from all its source blocks (removes itself
     * as an observer).
     */
    ~DependentBlock();

    /**
     * @brief Overwrites payload of an existing local parameter (no notifications).
     *
     * Calls Parameter::overwrite_payload_from() on the destination entry.
     * Does not notify parameter observers nor block observers.
     */
    void assign(const LhaID& key, std::shared_ptr<Parameter> param);

    /**
     * @brief Assigns a new value to an existing parameter (no notifications).
     *
     * Similar to @ref Block::assign(const LhaID&, scalar_t),
     * but does not trigger a cascade of updates.
     *
     * @param key   The LHA ID of the parameter.
     * @param value The new value to assign.
     */
    void assign(const LhaID& key, double value);

    /**
     * @brief Returns the source blocks of this dependent block.
     *
     * Overrides @ref Block::get_source_blocks().
     *
     * @return Map of source block names to their shared pointers.
     */
    std::unordered_map<std::string, std::shared_ptr<Block>> get_source_blocks() const override;

    /**
     * @brief Clears dependencies above this block.
     *
     * Unregisters this DependentBlock from all source blocks (removes it
     * from their observer lists).
     */
    void clear_above() override;

    /**
     * @brief Clears dependencies below this block (block graph only).
     *
     * Detaches from sources (clear_above()) and then calls clear_below() on
     * dependent blocks (observers). This function does not erase local parameters;
     * use destroy() for full teardown.
     */
    void clear_below() override;

    /**
     * @brief Destroys this DependentBlock and its dependents.
     *
     * Logs a debug message, clears links to source blocks, clears all local
     * parameters (after clearing their dependency trees), and recursively
     * destroys observer blocks.
     */
    void destroy() override;

    /**
     * @brief Marks the block content as invalid and propagates downstream.
     *
     * If already dirty, does nothing. Otherwise sets dirty=true and recursively
     * marks DependentBlock observers dirty as well.
     */
    void mark_dirty();

    /**
     * @brief Ensures this block is up-to-date (lazy recomputation).
     *
     * Calls the internal implementation that:
     *  - returns immediately if frozen or not dirty,
     *  - ensures upstream DependentBlock sources are up-to-date,
     *  - calls the recomputation lambda to rebuild/update local parameters,
     *  - then notifies parameter observers for each stored parameter.
     */
    void ensure_up_to_date() override { ensure_up_to_date_impl(); }

private:
    void ensure_up_to_date_impl();

    std::weak_ptr<DependentBlock> self;                                     ///< Self-reference for observer management.
    std::unordered_map<std::string, std::shared_ptr<Block>> sourceBlocks;   ///< Source blocks for dependencies.
    DepUpdateFunc recalculateLambda;                                        ///< Function used to recalculate this block's content.
    bool frozen = false;                                                    ///< Indicates if the block is frozen (no update).
    bool update_at_unfreeze = false;                                        ///< Indicates if an update is pending after unfreezing.
    bool dirty = true;                                                      ///< True if cached value is invalid and must be recomputed on demand.
};

#endif