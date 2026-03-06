#ifndef BLOCK_ACCESSOR_H
#define BLOCK_ACCESSOR_H

#include <functional>
#include <unordered_map>

#include "Include.h"
#include "Block.h"
#include "DependentParameter.h"

/**
 * @file BlockAccessor.h
 * @brief Provides an accessor for multiple parameter blocks.
 *
 * This file defines the BlockAccessor class, which manages multiple blocks
 * and provides a unified interface for accessing and modifying parameters
 * across blocks.
 */

/**
 * @class BlockAccessor
 * @brief A container / façade over multiple parameter blocks.
 *
 * This class derives from
 * `std::unordered_map<BlockName, std::shared_ptr<Block>>` and provides
 * convenience methods to:
 *
 *   - check whether a parameter exists in a given block,
 *   - get / set parameter values,
 *   - insert or update full Parameter objects,
 *   - query all values of a block,
 *   - remove parameters (with proper dependency cleanup through Block::remove),
 *   - query scales and sources (both at block and parameter level),
 *   - merge several BlockAccessor instances with or without priority.
 *
 * It is primarily used by higher-level parameter containers (e.g. a
 * Parameters class) to manipulate multiple physics blocks consistently.
 * @note Although the class publicly derives from an unordered_map, it also
 *       maintains an aliasing layer (alias_to_key_/key_to_name_) and stores
 *       blocks under a canonical key derived from block aliases.
 *       The public API should be used instead of raw map access.
 */
class BlockAccessor : public std::unordered_map<std::string, std::shared_ptr<Block>> {
public:

    using base_t = std::unordered_map<std::string, std::shared_ptr<Block>>;
    using base_t::operator[];

    /**
     * @brief Checks if a block exists with a given parameter.
     *
     * The function returns true if:
     *  - a block with name @p blockName exists in the accessor, and
     *  - this block contains a parameter with ID @p pdgCode.
     *
     * @param blockName The name of the block.
     * @param pdgCode   The LhaID (PDG-style identifier) of the parameter.
     * @return True if the parameter exists in the block, false otherwise.
     */
    bool has_param(const BlockName& blockName, LhaID pdgCode) const;

    /**
     * @brief Sets the value of a parameter in a specified block.
     *
     * Behavior:
     *  - If the block exists and already contains the parameter, its expected
     *    value is updated via Block::assign(id, value).
     *  - If the block exists but the parameter does not, a new Parameter is
     *    created with zero uncertainties and stored in the block.
     *  - If the block does not exist, an exception is thrown.
     *
     * @param blockName The name of the block.
     * @param pdgCode   The LhaID of the parameter.
     * @param value     The new central value to set.
     */
    void setValue(const BlockName& blockName, LhaID pdgCode, scalar_t value);

    // /**
    //  * @brief Sets the mode of a parameter in a specified block.
    //  * @param blockName The name of the block.
    //  * @param pdgCode   The PDG code of the parameter.
    //  * @param mode      The mode to set.
    //  */
    // void setMode(const std::string& blockName, LhaID pdgCode, ParameterMode mode);

    /**
     * @brief Retrieves the value of a parameter from a specified block.
     *
     * @param blockName The name of the block.
     * @param pdgCode   The LhaID of the parameter.
     * @return The current value of the parameter (including any shift).
     *
     * @throws std::invalid_argument if the block does not exist.
     *         If the block exists but the parameter does not, the underlying
     *         Block::retrieve may log and/or throw.
     */
    scalar_t getValue(const BlockName& blockName, LhaID pdgCode) const;

    /**
     * @brief Retrieves a parameter object from a specified block.
     *
     * @param blockName The name of the block.
     * @param id        The LhaID of the parameter.
     * @return Shared pointer to the Parameter.
     *
     * @throws std::invalid_argument if the block does not exist.
     *         If the block exists but the parameter does not, the underlying
     *         Block::retrieve may log and/or throw.
     */
    std::shared_ptr<Parameter> getParameter(const BlockName& blockName, LhaID id) const;

    /**
     * @brief Inserts or replaces a Parameter in a specified block.
     *
     * Calls Block::store_or_assign on the underlying block. If the block
     * does not exist, an exception is thrown.
     *
     * @param blockName The name of the block.
     * @param id        The LhaID of the parameter.
     * @param source    Shared pointer to the source Parameter.
     *
     * @throws std::invalid_argument if the block does not exist.
     */
    void setParameter(const BlockName& blockName, LhaID id, std::shared_ptr<Parameter> source);

    /**
     * @brief Retrieves all values from a specified block.
     *
     * Iterates over all parameters in the block and returns a map from
     * LhaID to current value (Parameter::get_val()).
     *
     * @param blockName The name of the block.
     * @return A map of LhaID to scalar_t values.
     *
     * @throws std::invalid_argument if the block does not exist.
     */
    std::map<LhaID, scalar_t> getAllValues(BlockName blockName);

    /**
     * @brief Retrieves the set of all block names present in the accessor.
     *
     * @return An unordered_set of BlockName keys.
     * @note This returns the set of BlockName objects tracked in key_to_name_
     *       (alias table). Direct mutation of the underlying unordered_map base
     *       may desynchronize these structures
     */
    std::unordered_set<BlockName> get_block_names() const;

    /**
     * @brief Removes a parameter from a specified block.
     *
     * Calls Block::remove(id) on the underlying block. If the block does not
     * exist, a warning is logged and no action is taken.
     *
     * @param block_name The name of the block.
     * @param id         The LhaID of the parameter to remove.*
     * @note Only the BlockName overload logs a warning and returns if the block
    *       does not exist. The string/string_view overloads call at(...) and throw.
     */
    void remove_item(const BlockName& block_name, LhaID id);

    /**
     * @brief Checks if a given block exists (by key or any alias).
     *
     * This hides the base-class contains(key) logic because BlockName has
     * custom comparison; it uses a helper to scan keys.
     *
     * @param block_name The block name to test.
     * @return True if the accessor contains a block with this name.
     */
    bool contains(const BlockName& block_name) const;

    /**
     * @brief Retrieves the scale associated with a given block.
     *
     * Internally forwards to Block::get_scale() after checking that:
     *  - the block exists, and
     *  - the block has a scale set (via Block::has_scale()).
     *
     * @param block_name The name of the block.
     * @return The scale value.
     *
     * @throws If the block does not exist or has no scale. Errors are logged.
     */
    double get_scale(const BlockName& block_name) const;

    /**
     * @brief Checks whether a specific block has an associated scale.
     *
     * @param block_name The name of the block.
     * @return True if the block exists and has a scale, false otherwise.
     *
     * @throws If the block does not exist (error is logged).
     */
    bool has_scale(const BlockName& block_name) const;

    /**
     * @brief Returns the source blocks of a given block.
     *
     * This is a thin wrapper around Block::get_source_blocks() for the
     * block named @p block_name.
     *
     * @param block_name The name of the block.
     * @return A map of source block names to their shared pointers.
     *
     * @throws If the block does not exist. Errors are logged.
     */
    std::unordered_map<std::string, std::shared_ptr<Block>> get_block_sources(const BlockName& block_name) const;

    /**
     * @brief Returns the source parameters contributing to a given parameter.
     *
     * This forwards to Parameter::get_source_parameters() for the parameter
     * identified by (block_name, id).
     *
     * @param block_name The name of the block.
     * @param id         The LhaID of the parameter.
     * @return A map of ParamId to shared_ptr<Parameter> describing the sources.
     *
     * @throws If the block or parameter does not exist. Errors are logged.
     */
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_parameter_sources(const BlockName& block_name, LhaID id) const;

    /**
     * @brief Recursively collects all "root" source parameters for a given set.
     *
     * For each input ParamId, this method:
     *  - follows parameter-level sources (Parameter::get_source_parameters()),
     *  - follows block-level sources (Block::get_source_blocks()),
     *  - recurses until no further sources are found.
     *
     * Any parameter ID that cannot be resolved in the current accessor or
     * has no sources is added to the result set.
     *
     * @param param_ids Set of parameter IDs from which to start the search.
     * @return Set of all leaf / root ParamId sources.
     */
    std::unordered_set<ParamId>
    get_all_source_parameters(const std::unordered_set<ParamId>& param_ids) const;
    
    /**
     * @brief Accessor to a block by name (alias-aware).
     *
     * This performs a linear search over the underlying map, using BlockName
     * equality semantics. If the block is not found, the current accessor
     * is printed and an error is logged before throwing.
     *
     * @param block_name Name of the block.
     * @return Reference to the shared_ptr<Block> stored under this name.
     */
    std::shared_ptr<Block>& at(const BlockName& block_name);

    /**
     * @brief Accessor to a block by name (const overload).
     *
     * Same semantics as the non-const overload.
     *
     * @param block_name Name of the block.
     * @return Const reference to the shared_ptr<Block> stored under this name.
     */
    const std::shared_ptr<Block>& at(const BlockName& block_name) const;

    /**
     * @brief Merges two BlockAccessor instances without resolving conflicts.
     *
     * For each block in @p lhs and @p rhs, a deep copy is made via the
     * `Block(std::shared_ptr<Block>)` copy constructor, and stored into
     * the resulting accessor.
     *
     * If both accessors contain a block with the same name, an error is
     * logged and the operation is considered invalid.
     *
     * @param lhs First BlockAccessor instance.
     * @param rhs Second BlockAccessor instance.
     * @return A shared pointer to a new BlockAccessor containing all unique blocks.
     *
     * @throws If overlapping blocks are found (error is logged and the call
     *         terminates via LOG_ERROR mechanics).
     */
    friend std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Merges two BlockAccessor instances with priority given to @p lhs.
     *
     * Semantics:
     *  - All blocks from @p rhs are copied into the result first.
     *  - Blocks from @p lhs are then merged in:
     *      * if a block exists only in lhs, it is copied in whole;
     *      * if a block exists in both, parameters from lhs are merged on top:
     *          - parameters only in lhs are added;
     *          - parameters in both are combined, with lhs central value kept,
     *            but uncertainties optionally taken from rhs if lhs has none.
     *
     * @param lhs First BlockAccessor instance (priority).
     * @param rhs Second BlockAccessor instance.
     * @return A shared pointer to a new BlockAccessor containing merged blocks.
     * @note Parameters present only in rhs are kept as-is (rhs is copied first).
     *       The merge loop iterates over lhs IDs only.
     */
    friend std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Extracts a subset of blocks (shallow copy of shared_ptrs).
     *
     * The returned accessor shares ownership of the same Block instances
     * (no deep copy is performed).
     *
     * @throws If a requested block does not exist (LOG_ERROR).
     */
    std::shared_ptr<BlockAccessor> operator[](std::unordered_set<BlockName> block_names);
    
    /**
     * @brief Stream output operator for BlockAccessor.
     *
     * Prints all block names and their parameter values using BlockAccessor::getAllValues:
     *
     * @code
     * Block SMINPUTS:
     *     id: value
     *     ...
     *
     * Block MASS:
     *     ...
     * @endcode
     *
     * @param os Output stream.
     * @param ba Shared pointer to the BlockAccessor to print.
     * @return The output stream.
     */
    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<BlockAccessor>);

    /**
     * @brief Inserts or replaces a block under a canonical key derived from its aliases.
     *
     * This updates alias tables (alias_to_key_, key_to_name_), stores the block
     * in the underlying map, and calls Block::bind_self() so the block can later
     * obtain a valid weak self pointer.
     *
     * @param name BlockName (with aliases).
     * @param blk  Block instance.
     */
    void emplace(const BlockName& name, std::shared_ptr<Block> blk);
    
    //TODO: docstring
    void detach_block(const BlockName& block_name);

    void detach_parameter(const BlockName& block_name, LhaID id);

    /**
     * @brief Erases a block by name/alias and removes all associated aliases.
     *
     * If the block cannot be resolved, this is a no-op.
     *
     * @param name Block name/alias.
     */
    void erase_block(const BlockName& name);

    /**
     * @name Convenience overloads (std::string / const char*)
     * @brief Thin wrappers around the std::string_view / BlockName API.
     *
     * These overloads exist purely for ergonomics so callers can pass
     * `std::string` or C-strings directly. They forward to the corresponding
     * `std::string_view` or `BlockName` implementations.
     *
     * @note All semantics (alias resolution, exceptions, logging) are defined by the
     *       forwarded-to overloads. In particular:
     *       - at(...) will throw std::invalid_argument if the block cannot be resolved;
     *       - contains(...) returns false if the name/alias cannot be resolved;
     *       - erase_block(...) is a no-op if the name/alias cannot be resolved.
     *
     * @warning The `std::string_view` overloads do not own storage. These wrappers
     *          are safe because they create the string_view from arguments whose
     *          lifetime extends through the call.
     */
    ///@{

    /// Erase by name/alias (std::string overload).
    void erase_block(const std::string& name) { return erase_block(std::string_view{name}); }
    /// Erase by name/alias (C-string overload).
    void erase_block(const char* name) { return erase_block(std::string_view{name}); }

    /// Alias-aware existence check (std::string overload).
    bool contains(const std::string& name) const { return contains(std::string_view{name}); }
    /// Alias-aware existence check (C-string overload).
    bool contains(const char* name) const { return contains(std::string_view{name}); }

    /// Alias-aware block access (std::string overload).
    std::shared_ptr<Block>& at(const std::string& name) { return at(std::string_view{name}); }
    /// Alias-aware block access (std::string overload, const).
    const std::shared_ptr<Block>& at(const std::string& name) const { return at(std::string_view{name}); }
    /// Alias-aware block access (C-string overload).
    std::shared_ptr<Block>& at(const char* name) { return at(std::string_view{name}); }
    /// Alias-aware block access (C-string overload, const).
    const std::shared_ptr<Block>& at(const char* name) const { return at(std::string_view{name}); }

    /// Insert/replace a block using a raw string name (converted to BlockName).
    void emplace(const std::string& name, std::shared_ptr<Block> blk) { emplace(BlockName(name), std::move(blk)); }
    /// Insert/replace a block using a C-string name (converted to BlockName).
    void emplace(const char* name, std::shared_ptr<Block> blk) { emplace(BlockName(name), std::move(blk)); }

    /// Parameter existence check (std::string overload).
    bool has_param(const std::string& block, const LhaID& id) const { return has_param(std::string_view{block}, id); }
    /// Parameter existence check (C-string overload).
    bool has_param(const char* block, const LhaID& id) const { return has_param(std::string_view{block}, id); }

    /// Value getter (std::string overload).
    scalar_t getValue(const std::string& block, const LhaID& id) const { return getValue(std::string_view{block}, id); }
    /// Value getter (C-string overload).
    scalar_t getValue(const char* block, const LhaID& id) const { return getValue(std::string_view{block}, id); }

    /// Value setter (std::string overload).
    void setValue(const std::string& block, const LhaID& id, scalar_t v) { setValue(std::string_view{block}, id, v); }
    /// Value setter (C-string overload).
    void setValue(const char* block, const LhaID& id, scalar_t v) { setValue(std::string_view{block}, id, v); }

    /// Remove parameter (std::string overload).
    void remove_item(const std::string& block, const LhaID& id) { remove_item(std::string_view{block}, id); }
    /// Remove parameter (C-string overload).
    void remove_item(const char* block, const LhaID& id) { remove_item(std::string_view{block}, id); }

    /// Scale existence check (std::string overload).
    bool has_scale(const std::string& block) const { return has_scale(std::string_view{block}); }
    /// Scale existence check (C-string overload).
    bool has_scale(const char* block) const { return has_scale(std::string_view{block}); }

    /// Scale getter (std::string overload).
    double get_scale(const std::string& block) const { return get_scale(std::string_view{block}); }
    /// Scale getter (C-string overload).
    double get_scale(const char* block) const { return get_scale(std::string_view{block}); }

    ///@}
    
private:
    
    bool has_scale(std::string_view block) const;
    double get_scale(std::string_view block) const;
    
    bool contains(std::string_view name) const;
    void emplace(std::string_view name, std::shared_ptr<Block> blk) {
        emplace(BlockName(std::string(name)), std::move(blk));
    }
    std::shared_ptr<Block>& at(std::string_view name);
    const std::shared_ptr<Block>& at(std::string_view name) const;
    
    bool has_param(std::string_view block, const LhaID& id) const;
    scalar_t getValue(std::string_view block, const LhaID& id) const;
    void setValue(std::string_view block, const LhaID& id, scalar_t value);
    void remove_item(std::string_view block, const LhaID& id);
    void erase_block(std::string_view alias_or_key);
    
    static std::unordered_set<std::string> norm_aliases(const BlockName& n);
    static std::string choose_key(const std::unordered_set<std::string>& aliases_norm);
    static std::string normalize(std::string_view s);
    
    std::string key_for(std::string_view alias) const;
    std::string key_for(const BlockName& name) const;
    std::string key_for(const std::string& name) const;
    
    std::string resolve_key(std::string_view name) const;
    std::string resolve_key(const BlockName& name) const;
    void merge_name_into_key(const std::string& key, const BlockName& name);
    
    /// Maps normalized alias -> canonical key used in the underlying unordered_map.
    std::unordered_map<std::string, std::string> alias_to_key_;

    /// Maps canonical key -> full BlockName (with all known aliases merged).
    std::unordered_map<std::string, BlockName> key_to_name_;
};

#endif