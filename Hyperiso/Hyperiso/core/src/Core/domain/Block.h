#ifndef BLOCK_H
#define BLOCK_H

#include <array>
#include <functional>

#include "Include.h"
#include "Parameter.h"
#include "IStorage.h"

/**
 * @file Block.h
 * @brief Defines classes used to store parameters and to build derived/dependent parameter blocks.
 *
 * This file declares two main abstractions:
 * - @ref Block: a generic container of @ref Parameter objects indexed by @ref LhaID,
 * - @ref DependentBlock: a derived block whose content is recomputed from one or more source blocks.
 *
 * These classes form one of the core layers of the dependency/update system:
 * - a block stores parameters,
 * - parameters may notify dependent parameters,
 * - blocks may notify dependent blocks,
 * - dependent blocks may lazily recompute their content from upstream blocks.
 *
 * Main features provided here:
 * - storage/retrieval of parameters,
 * - optional block-wide scale,
 * - observer registration between blocks,
 * - block/parameter freeze-unfreeze semantics,
 * - explicit dependency clearing and destruction,
 * - lazy recomputation for dependent blocks.
 *
 * @see Parameter
 * @see DependentParameter
 * @see IStorage
 */

class BlockSrc;
class Block;
class DependentBlock;

/**
 * @typedef DepUpdateFunc
 * @brief Function type used to recompute a @ref DependentBlock from its sources.
 *
 * The function receives:
 * - a @ref BlockSrc view exposing the source blocks and their values,
 * - a shared pointer to the dependent block being recomputed.
 *
 * The callback is expected to update the dependent block in place, typically through
 * calls such as `store_or_assign(...)`.
 *
 * @see BlockSrc
 * @see DependentBlock
 */
typedef std::function<void(const BlockSrc&, std::shared_ptr<DependentBlock>)> DepUpdateFunc;


/**
 * @typedef DepUpdateFunc
 * @brief Function type used to recompute a @ref DependentBlock from its sources.
 *
 * The function receives:
 * - a @ref BlockSrc view exposing the source blocks and their values,
 * - a shared pointer to the dependent block being recomputed.
 *
 * The callback is expected to update the dependent block in place, typically through
 * calls such as `store_or_assign(...)`.
 *
 * @see BlockSrc
 * @see DependentBlock
 */
class Block : public IStorage<LhaID, Parameter>, public std::enable_shared_from_this<Block> {
public:
    /// Name of the block (e.g. "SMINPUTS", "MASS", "FWCOEF"...)
    BlockName blockname {""};

    /**
     * @brief Default constructor.
     *
     * Constructs an empty block with no parameters, no observers, and no scale.
     */
    Block() = default;

    /**
     * @brief Copy-like constructor from another block.
     *
     * This constructor copies the content and metadata of @p other by delegating
     * to @ref copy(). Stored parameter pointers are copied as shared pointers.
     *
     * After the copy, parameter owner-block links are rebound to this new block
     * whenever a valid self reference is available.
     *
     * @param other Block to copy from.
     */
    Block(std::shared_ptr<Block> other);

    /**
     * @brief Returns the name of the block.
     * @return Block name.
     */
    BlockName get_name() const { return this->blockname; }

    /**
     * @brief Stores a new parameter under the given LHA id.
     *
     * If the key already exists, the parameter is not overwritten and a debug
     * message is emitted.
     *
     * If the block has been bound to itself through @ref bind_self(), the stored
     * parameter receives this block as owner block.
     *
     * @param id    LHA identifier of the parameter.
     * @param param Parameter object to store.
     */
    void store(const LhaID& id, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Replaces the payload of an existing parameter from another parameter object.
     *
     * This does not replace the shared pointer stored in the block. Instead, it copies
     * the numerical/configuration payload into the already stored destination parameter
     * using `Parameter::overwrite_payload_from(...)`.
     *
     * Current notification behavior:
     * - the destination parameter notifies its own observers,
     * - block-level observers are reached indirectly through the parameter/owner-block chain,
     * - no extra explicit `Block::notifyObservers()` call is performed here.
     *
     * @param key   Identifier of the parameter to overwrite.
     * @param param Source parameter whose payload is copied.
     */
    void assign(const LhaID& key, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Sets the expected value of an existing parameter.
     *
     * This calls `Parameter::set_expected(value)` on the stored parameter.
     *
     * Current notification behavior:
     * - notifications are delegated to the parameter layer,
     * - no extra explicit `Block::notifyObservers()` is triggered here.
     *
     * @param key   Identifier of the parameter to update.
     * @param value New value.
     */
    void assign(const LhaID& key, scalar_t value);

    /**
     * @brief Stores a parameter if absent, otherwise assigns into the existing one.
     *
     * This is a convenience method:
     * - calls @ref store() if the key is not present,
     * - calls @ref assign(const LhaID&, std::shared_ptr<Parameter>) otherwise.
     *
     * @param key   Identifier of the parameter.
     * @param param Parameter to store or assign.
     */
    void store_or_assign(const LhaID& key, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Checks whether the block contains a parameter with the given id.
     *
     * Before the lookup, the block is forced up-to-date through
     * `ensure_up_to_date()`. For a plain @ref Block this is a no-op, but for
     * @ref DependentBlock this may trigger lazy recomputation.
     *
     * @param key LHA identifier.
     * @return True if the parameter exists.
     */
    bool contains(const LhaID& key) const override;

    /**
     * @brief Retrieves a stored parameter by id.
     *
     * The block is first ensured up-to-date.
     * If the id is missing, an error is logged.
     *
     * @param id LHA identifier to retrieve.
     * @return Shared pointer to the stored parameter.
     */
    std::shared_ptr<Parameter> retrieve(const LhaID& id) override;

    /**
     * @brief Removes a parameter and clears its downstream dependency subtree.
     *
     * The parameter is first asked to clear its own downstream dependencies via
     * `Parameter::clear_below()`, then removed from the block.
     *
     * After that, all observer blocks currently attached to this block are destroyed
     * recursively.
     *
     * @param key Identifier of the parameter to remove.
     */
    void remove(const LhaID& key) override;

    /**
     * @brief Returns the set of all parameter ids stored in the block.
     *
     * The block is ensured up-to-date before extracting the ids.
     *
     * @return Set of stored ids.
     */
    std::unordered_set<LhaID> getAllIDs();

    /**
     * @brief Returns the internal parameter map.
     *
     * The block is ensured up-to-date before returning the map.
     *
     * @return Const reference to the internal storage map.
     */
    const std::map<LhaID, std::shared_ptr<Parameter>>& getItems();

    /**
     * @brief Sets the owner ParameterType on all contained parameters.
     *
     * This is forwarded to each stored parameter via `Parameter::set_owner(type)`.
     *
     * @param type Owner type to apply.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Adds a block observer if not already present.
     *
     * Observer uniqueness is checked by pointer identity.
     *
     * @param observer Observer block to add.
     */
    void addObserver(std::shared_ptr<Block> observer);

    /**
     * @brief Removes a previously registered observer block.
     *
     * Removal is performed by pointer identity.
     *
     * @param observer Observer block to remove.
     */
    void removeObserver(std::shared_ptr<Block> observer);

    /**
     * @brief Notifies all observer blocks that this block has changed.
     *
     * For each non-null observer, `observer->update()` is called.
     * Null observer entries are cleaned afterward.
     */
    void notifyObservers();

    /**
     * @brief Returns the current list of observer blocks.
     * @return Observer list.
     */
    std::vector<std::shared_ptr<Block>> getObservers() const;

    /**
     * @brief Initialization hook.
     *
     * Plain blocks do nothing by default. Derived classes may override this
     * to register dependencies or pre-fill content.
     */
    virtual void init() {};
    
    /**
     * @brief Updates the block.
     *
     * Default behavior: calls `update()` on every contained parameter.
     *
     * For plain blocks this is an eager propagation method. Dependent blocks
     * override this with lazy dirty-marking semantics.
     */
    virtual void update();

    /**
     * @brief Destroys the block and its downstream block dependencies.
     *
     * Current behavior:
     * - snapshots and detaches block observers,
     * - clears parameter-side downstream dependencies,
     * - clears local storage,
     * - recursively destroys dependent observer blocks.
     *
     * @note This is explicit graph teardown logic; it is not automatically called
     * by the destructor.
     */
    virtual void destroy();
    
    /**
     * @brief Freezes all parameters contained in the block.
     *
     * This delegates to `Parameter::freeze()` on each stored parameter.
     */
    virtual void freeze();

    /**
     * @brief Unfreezes all parameters contained in the block.
     *
     * This delegates to `Parameter::unfreeze()` on each stored parameter.
     */
    virtual void unfreeze();

    /**
     * @brief Copies the content and metadata from another block.
     *
     * Copied elements:
     * - storage map,
     * - block name,
     * - scale if present.
     *
     * Stored parameter pointers are copied as shared pointers. If this block is already
     * self-bound, copied parameters have their owner-block link rebound to this block.
     *
     * @param other Source block.
     */
    void copy(std::shared_ptr<Block> other);

    /**
     * @brief Deep-copies this block as a plain independent block.
     *
     * All contained Parameter payloads are duplicated, observer links and
     * owner-block links are not shared with the source block. Dependency graph
     * objects are intentionally not cloned here; dependent blocks are rebuilt
     * by the owning Parameters strategy inside each runtime context.
     */
    std::shared_ptr<Block> deep_clone_plain() const;

    /**
     * @brief Clears parameter dependencies above this block.
     *
     * Delegates `clear_above()` to all contained parameters.
     */
    virtual void clear_above();

    /**
     * @brief Clears parameter dependencies below this block.
     *
     * A snapshot of current parameter pointers is built first to avoid iterator
     * invalidation during recursive dependency clearing.
     */
    virtual void clear_below();

    /**
     * @brief Returns whether the block has an associated scale.
     * @return True if a scale is set.
     */
    bool has_scale();

    /**
     * @brief Sets the block-wide scale.
     *
     * @param scale Scale to assign.
     *
     * @note Setting the scale twice is considered a logic error.
     */
    void set_scale(double scale);

    /**
     * @brief Returns the block-wide scale.
     * @return The scale value.
     *
     * @note Accessing an unset scale is considered a logic error.
     */
    double get_scale();

    /**
     * @brief Erases an entry from the local storage map only.
     *
     * Unlike @ref remove(), this does not clear any dependency relationships and
     * does not destroy downstream blocks. It is a pure local erase.
     *
     * @param id Identifier to erase.
     */
    void erase_local(const LhaID& id);

    /**
     * @brief Returns the source blocks of this block.
     *
     * Plain blocks have no source blocks, so the default implementation returns
     * an empty map.
     *
     * @return Source block map.
     */
    virtual std::unordered_map<std::string, std::shared_ptr<Block>> get_source_blocks() const;

    /**
     * @brief Ensures that the block content is up to date.
     *
     * For a plain @ref Block this is a no-op.
     * @ref DependentBlock overrides it to implement lazy recomputation.
     */ 
    virtual void ensure_up_to_date() {}

    /**
     * @brief Binds an explicit self shared_ptr to the block.
     *
     * This is used to make owner-block links available even in situations where
     * `shared_from_this()` is not yet safe or convenient.
     *
     * After binding, all currently stored parameters are rebound to this block
     * as their owner block.
     *
     * @param self Shared pointer referring to this block instance.
     */
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
     * The destructor itself does not perform dependency teardown.
     * Use @ref destroy() when graph-level destruction is required.
     */
    ~Block() {}

    /**
     * @brief Stream output operator for a block.
     *
     * Prints block content in a human-readable form:
     * @code
     * Block MASS:
     *     6: 172.5
     *     24: 80.379
     * @endcode
     *
     * @param os Output stream.
     * @param ba Block to print.
     * @return Output stream.
     */
    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<Block>);
protected:
    std::vector<std::shared_ptr<Block>> observers;      /// List of observing blocks notified through @ref notifyObservers().
    std::map<LhaID, std::shared_ptr<Parameter>> items;  /// Internal storage of parameters indexed by LHA id.
    std::optional<double> scale;                        /// Optional block-wide scale.

    std::weak_ptr<Block> self_;                         /// Optional explicit self-reference used to rebind owner_block on contained parameters.

    /**
     * @brief Returns the best available weak self-reference.
     *
     * Priority:
     * - explicit binding through @ref bind_self(),
     * - fallback to `shared_from_this()` when available,
     * - empty weak_ptr otherwise.
     *
     * @return Weak self-reference.
     */
    std::weak_ptr<Block> self_weak() {
        if (auto s = self_.lock()) return self_;
        try { return shared_from_this(); } catch (...) { return {}; }
    }
};

/**
 * @class DependentBlock
 * @brief Block whose content is derived from one or more source blocks.
 *
 * A DependentBlock is a lazily recomputed block:
 * - it observes one or more source blocks,
 * - when sources change, it becomes dirty,
 * - actual recomputation is deferred until a read-style access forces it
 *   (e.g. @ref retrieve(), @ref contains(), @ref getItems(), @ref getAllIDs()).
 *
 * The actual recomputation logic is provided by a user callback of type
 * @ref DepUpdateFunc.
 *
 * Additional features:
 * - freeze/unfreeze support with deferred update,
 * - detach/reattach support to temporarily disable the dependency graph,
 * - recursive dirty propagation to downstream dependent blocks and parameters.
 *
 * @see Block
 * @see DepUpdateFunc
 * @see BlockSrc
 */
class DependentBlock : public Block {
public:
    /**
     * @brief Constructs a dependent block from source blocks and a recomputation callback.
     *
     * The constructor only stores the dependency data. Actual observer registration
     * is performed later by @ref init().
     *
     * @param sources Source blocks used as inputs.
     * @param recalculateFunc Callback used to rebuild/update this block.
     */
    explicit DependentBlock(
    const std::unordered_map<std::string, std::shared_ptr<Block>>& sources,
    DepUpdateFunc recalculateFunc)
    : sourceBlocks(sources),
      recalculateLambda(recalculateFunc),
      saved_sourceBlocks(sources),
      saved_recalculateLambda(recalculateFunc),
      frozen(false) {}

    /**
     * @brief Checks whether this dependent block uses a given source block name.
     *
     * @param blockName Source block name.
     * @return True if the source exists and is non-null.
     */
    bool dependsOn(const std::string& blockName);

    /**
     * @brief Initializes the dependency graph for this block.
     *
     * Registers this block as observer of all current source blocks.
     * Must be called before the dependency graph is expected to function.
     */
    void init() override;

    /**
     * @brief Marks the block as dirty.
     *
     * Current behavior:
     * - if frozen, records that an update must happen upon unfreeze,
     * - otherwise propagates dirtiness downstream using @ref mark_dirty().
     *
     * This method does not recompute immediately.
     */
    void update() override;

    /**
     * @brief Freezes the block and all contained parameters.
     *
     * While frozen, updates are deferred.
     */
    void freeze() override;

    /**
     * @brief Unfreezes the block and contained parameters.
     *
     * If an update was requested while frozen, the block is marked dirty again
     * upon unfreeze.
     */
    void unfreeze() override;

    /**
     * @brief Temporarily detaches this block from its dependency graph.
     *
     * Current behavior:
     * - ensures the block is up to date before detaching,
     * - unregisters this block from its source observers,
     * - detaches dependent parameters when applicable,
     * - clears source/recalculation state from the active dependency graph,
     * - preserves a saved copy of the original dependency definition,
     * - notifies downstream observers.
     *
     * After detachment, the block behaves like a frozen snapshot of its current content.
     *
     * @see reattach
     */
    void detach();

    /**
     * @brief Reattaches a previously detached dependent block.
     *
     * Restores:
     * - source blocks,
     * - recomputation callback,
     * - parameter dependency state,
     * - observer registration to sources.
     *
     * The block is then marked dirty again so that future accesses recompute
     * from restored sources.
     *
     * @see detach
     */
    void reattach();

    /**
     * @brief Destructor.
     *
     * Unregisters this dependent block from its source blocks when still attached.
     */
    ~DependentBlock();

    /**
     * @brief Overwrites the payload of a local parameter without triggering block-level notification.
     *
     * This is a lighter-weight internal assign than @ref Block::assign(...):
     * the destination parameter payload is overwritten, but this method does not
     * explicitly notify downstream blocks.
     *
     * @param key   Identifier to overwrite.
     * @param param Source parameter payload.
     */
    void assign(const LhaID& key, std::shared_ptr<Parameter> param);

    /**
     * @brief Assigns a new scalar value to a local parameter.
     *
     * This is the dependent-block-local variant used during recomputation.
     *
     * @param key   Identifier to update.
     * @param value New value.
     */
    void assign(const LhaID& key, double value);

    /**
     * @brief Returns the map of source blocks for this dependent block.
     * @return Source block map.
     */
    std::unordered_map<std::string, std::shared_ptr<Block>> get_source_blocks() const override;

    /**
     * @brief Clears upstream block dependencies.
     *
     * Unregisters this block from all current source blocks.
     */
    void clear_above() override;

    /**
     * @brief Clears downstream block dependencies.
     *
     * This first detaches from sources, then recursively asks observing blocks to
     * clear their own downstream dependencies.
     */
    void clear_below() override;

    /**
     * @brief Destroys this dependent block and its downstream block dependencies.
     *
     * Current behavior:
     * - detaches from sources,
     * - clears downstream parameter dependencies,
     * - clears local items,
     * - recursively destroys observing blocks.
     */
    void destroy() override;

    /**
     * @brief Marks this block dirty and propagates the dirty state downstream.
     *
     * If the block is already dirty, nothing is done.
     *
     * Propagation rules:
     * - downstream @ref DependentBlock observers receive `mark_dirty()`,
     * - downstream dependent parameters are notified through
     *   `Parameter::notifyParamObserversOnly()`.
     */
    void mark_dirty();

    /**
     * @brief Ensures this dependent block is up to date.
     *
     * If the block is dirty and not frozen, upstream dependent blocks are first
     * ensured up to date, then the recomputation callback is executed, and finally
     * downstream parameter/block observers are notified in a controlled way.
     */
    void ensure_up_to_date() override { ensure_up_to_date_impl(); }

private:
    /**
     * @brief Internal lazy recomputation routine.
     *
     * This method implements the actual on-demand update policy for dependent blocks.
     */
    void ensure_up_to_date_impl();

    std::weak_ptr<DependentBlock> self;                                     /// Weak self-reference specialized as DependentBlock for source observer management.
    std::unordered_map<std::string, std::shared_ptr<Block>> sourceBlocks;   /// Active source blocks currently used by the dependency graph.
    DepUpdateFunc recalculateLambda;                                        /// Active recomputation callback.
    bool frozen = false;                                                    ///< Indicates if the block is frozen (no update).
    bool update_at_unfreeze = false;                                        /// True if an update request happened while frozen and should be replayed after unfreeze.
    bool dirty = true;                                                      /// True if the cached content is invalid and must be recomputed on next access.

    /**
     * @brief Saved source blocks used to restore dependencies after @ref detach().
     *
     * These are not necessarily the currently active source blocks. They represent
     * the dependency definition that should be restored by @ref reattach().
     */
    std::unordered_map<std::string, std::shared_ptr<Block>> saved_sourceBlocks;

    /**
     * @brief Saved recomputation callback used to restore dependencies after @ref detach().
     *
     * This mirrors @ref saved_sourceBlocks and allows the dependent block to resume
     * its original update logic after reattachment.
     */
    DepUpdateFunc saved_recalculateLambda;

    /**
     * @brief Indicates whether the dependency graph is currently detached.
     *
     * When true:
     * - the block no longer observes upstream blocks,
     * - active sources/callback may be cleared,
     * - the block behaves as a detached snapshot until @ref reattach() is called.
     */
    bool dependency_detached = false;
};

#endif