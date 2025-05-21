/**
 * @file Block.h
 * @brief Defines classes for managing parameter blocks and their dependencies.
 *
 * This file declares a hierarchy of classes used to store, manage, and update parameter values across different
 * parameter instances. It includes support for block dependencies and automatic notification of updates.
 */
#if !defined(BLOCK_H)
#define BLOCK_H

#include <array>
#include <functional>

#include "Include.h"
#include "Parameter.h"
#include "IStorage.h"

class Block;
class DependentBlock;

/**
 * @typedef DepUpdateFunc
 * @brief Defines the function signature used to update dependent blocks.
 *
 * The function takes a map of source blocks and a shared pointer to the dependent block itself.
 */
typedef std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> DepUpdateFunc;

/**
 * @class Block
 * @brief A class for storing and managing a collection of parameters identified by LHA IDs.
 *
 * Supports storing, updating, freezing, and notifying dependent observers when parameters change.
 */
class Block : public IStorage<LhaID, Parameter> {
public:
    BlockName blockname {""}; ///< Name of the block.

    Block() = default;
    Block(std::shared_ptr<Block> other);

    BlockName get_name() const { return this->blockname; }

    /**
     * @brief Stores a new parameter.
     * @param id The LHA ID of the parameter.
     * @param param Shared pointer to the Parameter to store.
     */
    void store(const LhaID& id, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Assigns a new value to an existing parameter.
     * @param key The LHA ID of the parameter.
     * @param param Shared pointer to the new Parameter.
     */
    void assign(const LhaID& key, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Assigns a new value (double) to an existing parameter.
     * @param key The LHA ID of the parameter.
     * @param value The new value to assign.
     */
    void assign(const LhaID& key, scalar_t value);

    /**
     * @brief Stores or assigns a parameter depending on whether it already exists.
     * @param key The LHA ID of the parameter.
     * @param param Shared pointer to the Parameter.
     */
    void store_or_assign(const LhaID& key, std::shared_ptr<Parameter> param) override;

    /**
     * @brief Checks if a parameter exists.
     * @param key The LHA ID of the parameter.
     * @return True if the parameter exists, false otherwise.
     */
    bool contains(const LhaID& key) const override;

    /**
     * @brief Retrieves a parameter.
     * @param id The LHA ID of the parameter.
     * @return Shared pointer to the Parameter.
     */
    std::shared_ptr<Parameter> retrieve(const LhaID& id) override;

    /**
     * @brief Removes a parameter.
     * @param key The LHA ID of the parameter to remove.
     */
    void remove(const LhaID& key) override;

    /**
     * @brief Retrieves all parameter IDs in the block.
     * @return Set of LHA IDs.
     */
    std::unordered_set<LhaID> getAllIDs();

    /**
     * @brief Retrieves all stored parameters.
     * @return Const reference to the internal map of items.
     */
    const std::map<LhaID, std::shared_ptr<Parameter>>& getItems() { return this->items; };

    /**
     * @brief Sets the owner type for all parameters in the block.
     * @param type The ParameterType to set as owner.
     */
    void set_owner(ParameterType type);

    /**
     * @brief Adds an observer block that will be notified on updates.
     * @param observer Shared pointer to the observer block.
     */
    void addObserver(std::shared_ptr<Block> observer);

    /**
     * @brief Removes an observer block.
     * @param observer Shared pointer to the observer block.
     */
    void removeObserver(std::shared_ptr<Block> observer);

    /**
     * @brief Notifies all observer blocks of an update.
     */
    void notifyObservers();

    /**
     * @brief Virtual method to update the block (default behavior: update all parameters).
     */
    virtual void update();

    /**
     * @brief Virtual method to freeze the block (prevent updates).
     */
    virtual void freeze();

    /**
     * @brief Virtual method to unfreeze the block (allow updates).
     */
    virtual void unfreeze();

    /**
     * @brief Copies the contents of another block into this one.
     * @param other Shared pointer to the other block.
     */
    void copy(std::shared_ptr<Block> other);

    /**
     * @brief Destructor. Notifies observers when the block is destroyed.
     */
    ~Block() { notifyObservers(); }

    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<Block>);
protected:
    std::vector<std::shared_ptr<Block>> observers;      ///< List of observer blocks.
    std::map<LhaID, std::shared_ptr<Parameter>> items;  ///< Map of parameters.
};

/**
 * @class DependentBlock
 * @brief A specialized Block that depends on other Blocks and updates automatically when they change.
 *
 * Implements dependency management and lazy update mechanisms (freeze/unfreeze behavior).
 */
class DependentBlock : public Block, public std::enable_shared_from_this<DependentBlock> {
public:
    /**
     * @brief Constructs a DependentBlock.
     * @param sources Map of source block names to block pointers.
     * @param recalculateFunc Function to recalculate the dependent block.
     */
    explicit DependentBlock(const std::unordered_map<std::string, std::shared_ptr<Block>>& sources, DepUpdateFunc recalculateFunc) 
        : sourceBlocks(std::move(sources)), recalculateLambda(std::move(recalculateFunc)), frozen(false) {}

    /**
     * @brief Checks if this block depends on a specific source block.
     * @param blockName Name of the source block.
     * @return True if dependency exists, false otherwise.
     */
    bool dependsOn(const std::string& blockName);

    /**
     * @brief Initializes dependency tracking (must be called after construction).
     */
    void init();

    /**
     * @brief Updates the dependent block.
     */
    void update() override;

    /**
     * @brief Freezes updates to the block.
     */
    void freeze() override;

    /**
     * @brief Unfreezes updates to the block and triggers an update if needed.
     */
    void unfreeze() override;

    /**
     * @brief Destructor. Cleans up dependencies.
     */
    ~DependentBlock();

    /**
     * @brief Assigns a new parameter to an existing entry.
     * @param key The LHA ID of the parameter.
     * @param param Shared pointer to the new Parameter.
     */
    void assign(const LhaID& key, std::shared_ptr<Parameter> param);

    /**
     * @brief Assigns a new value to an existing parameter.
     * @param key The LHA ID of the parameter.
     * @param value The new value to assign.
     */
    void assign(const LhaID& key, double value);

private:
    std::shared_ptr<DependentBlock> self;                                   ///< Self-reference for observer management.
    std::unordered_map<std::string, std::shared_ptr<Block>> sourceBlocks;   ///< Source blocks for dependencies.
    DepUpdateFunc recalculateLambda;                                        ///< Function used to recalculate this block's content.
    bool frozen;                                                            ///< Indicates if the block is frozen (no update).
    bool update_at_unfreeze;                                                ///< Indicates if an update is pending after unfreezing.
};

#endif