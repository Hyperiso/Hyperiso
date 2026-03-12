#ifndef BLOCK_ACCESSOR_H
#define BLOCK_ACCESSOR_H

#include <functional>
#include <unordered_map>

#include "Include.h"
#include "Block.h"
#include "DependentParameter.h"

/**
 * @file BlockAccessor.h
 * @brief Alias-aware façade for accessing and manipulating multiple parameter blocks.
 *
 * This file defines @ref BlockAccessor, a high-level container built on top of
 * a map of block names to @ref Block instances.
 *
 * Compared to using raw blocks directly, BlockAccessor adds:
 * - alias-aware block lookup,
 * - convenience getters/setters for parameter values and Parameter objects,
 * - support for dependency introspection (block sources and parameter sources),
 * - utilities to detach / reattach dependent blocks and dependent parameters,
 * - merge operators for combining several block collections,
 * - extraction of sub-accessors on a subset of blocks.
 *
 * The class publicly derives from
 * `std::unordered_map<std::string, std::shared_ptr<Block>>`, but callers are
 * expected to use the provided API instead of raw map operations because the
 * class maintains an internal alias-resolution layer:
 * - `alias_to_key_` maps normalized aliases to a canonical storage key,
 * - `key_to_name_` maps canonical keys back to the full @ref BlockName object.
 *
 * @see Block
 * @see DependentBlock
 * @see Parameter
 * @see DependentParameter
 */

/**
 * @class BlockAccessor
 * @brief Alias-aware container / façade over several parameter blocks.
 *
 * A BlockAccessor stores blocks under canonical internal keys while exposing a
 * more user-friendly interface based on @ref BlockName aliases.
 *
 * Main features:
 * - check whether a block or parameter exists,
 * - read/write scalar values,
 * - read/write full @ref Parameter objects,
 * - retrieve all values from a block,
 * - remove parameters with proper dependency cleanup,
 * - inspect source blocks / source parameters,
 * - recursively find leaf source parameters,
 * - merge several accessors,
 * - detach/reattach dependent blocks or dependent parameters.
 *
 * @note Although this class publicly inherits from `std::unordered_map`, direct
 *       manipulation of the base map may desynchronize alias metadata
 *       (`alias_to_key_`, `key_to_name_`). Prefer the dedicated API.
 */
class BlockAccessor : public std::unordered_map<std::string, std::shared_ptr<Block>> {
public:

    /// Underlying associative container type.
    using base_t = std::unordered_map<std::string, std::shared_ptr<Block>>;

    /// Re-expose the base `operator[]` for internal / advanced use.
    using base_t::operator[];

    /**
     * @brief Checks whether a given parameter exists inside a named block.
     *
     * This method first resolves the block name through the alias system, then
     * delegates to `Block::contains(id)`.
     *
     * @param blockName Name (or alias set) of the block.
     * @param pdgCode   LHA identifier of the parameter inside the block.
     * @return True if the block exists and contains the parameter, false otherwise.
     */

    bool has_param(const BlockName& blockName, LhaID pdgCode) const;

    /**
     * @brief Sets the value of a parameter inside an existing block.
     *
     * Behavior:
     * - if the block exists and the parameter already exists, delegates to
     *   `Block::assign(id, value)`,
     * - if the block exists but the parameter does not, creates a new
     *   @ref Parameter with zero uncertainties and stores it,
     * - if the block does not exist, throws `std::invalid_argument`.
     *
     * @param blockName Name of the target block.
     * @param pdgCode   LHA identifier of the parameter.
     * @param value     New scalar value.
     */
    void setValue(const BlockName& blockName, LhaID pdgCode, scalar_t value);

    /**
     * @brief Retrieves the current value of a parameter from a block.
     *
     * This is equivalent to:
     * @code
     * at(blockName)->retrieve(pdgCode)->get_val();
     * @endcode
     *
     * @param blockName Name of the block.
     * @param pdgCode   LHA identifier of the parameter.
     * @return Current parameter value.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     *         If the block exists but the parameter does not, the underlying
     *         `Block::retrieve()` may log and/or throw.
     */
    scalar_t getValue(const BlockName& blockName, LhaID pdgCode) const;

    /**
     * @brief Retrieves the full @ref Parameter object stored in a block.
     *
     * @param blockName Name of the block.
     * @param id        LHA identifier of the parameter.
     * @return Shared pointer to the stored parameter.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     *         If the block exists but the parameter does not, the underlying
     *         `Block::retrieve()` may log and/or throw.
     */
    std::shared_ptr<Parameter> getParameter(const BlockName& blockName, LhaID id) const;

    /**
     * @brief Inserts or updates a full @ref Parameter object in a block.
     *
     * Delegates to `Block::store_or_assign(id, source)` on the resolved block.
     *
     * @param blockName Name of the block.
     * @param id        LHA identifier of the parameter.
     * @param source    Shared pointer to the source parameter object.
     *
     * @throws std::invalid_argument if the block does not exist.
     */
    void setParameter(const BlockName& blockName, LhaID id, std::shared_ptr<Parameter> source);

    /**
     * @brief Returns all current scalar values stored in a block.
     *
     * Each stored parameter is converted to its current scalar value through
     * `Parameter::get_val()`.
     *
     * @param blockName Name of the block.
     * @return Ordered map of `LhaID -> scalar_t`.
     *
     * @throws std::invalid_argument if the block does not exist.
     */
    std::map<LhaID, scalar_t> getAllValues(BlockName blockName);

    /**
     * @brief Returns the set of all known block names.
     *
     * This returns the values stored in the alias metadata table
     * `key_to_name_`, not the raw keys of the underlying unordered_map.
     *
     * @return Set of block names with aliases preserved.
     */
    std::unordered_set<BlockName> get_block_names() const;

    /**
     * @brief Removes one parameter from a block.
     *
     * Delegates to `Block::remove(id)` which performs proper dependency cleanup.
     *
     * Behavior of this overload:
     * - if the block exists, the parameter is removed,
     * - if the block does not exist, a warning is logged and nothing happens.
     *
     * @param block_name Name of the block.
     * @param id         LHA identifier of the parameter to remove.
     */
    void remove_item(const BlockName& block_name, LhaID id);

    /**
     * @brief Checks whether a block exists (alias-aware).
     *
     * This method resolves aliases through the internal alias map instead of
     * using the raw underlying map lookup.
     *
     * @param block_name Name (or aliases) of the block.
     * @return True if the block exists, false otherwise.
     */
    bool contains(const BlockName& block_name) const;

    /**
     * @brief Returns the scale attached to a given block.
     *
     * Delegates to `Block::get_scale()` after resolving the block.
     *
     * @param block_name Name of the block.
     * @return Block scale.
     *
     * @throws If the block does not exist or if the block has no scale set.
     */
    double get_scale(const BlockName& block_name) const;

    /**
     * @brief Checks whether a given block has a scale attached.
     *
     * @param block_name Name of the block.
     * @return True if the block exists and has a scale, false otherwise.
     *
     * @throws If the block does not exist.
     */
    bool has_scale(const BlockName& block_name) const;

    /**
     * @brief Returns the source blocks of a given block.
     *
     * For plain @ref Block objects, this is usually an empty map.
     * For @ref DependentBlock objects, this returns the actual dependency map.
     *
     * @param block_name Name of the block.
     * @return Map of source block names to source block pointers.
     *
     * @throws If the block does not exist.
     */
    std::unordered_map<std::string, std::shared_ptr<Block>> get_block_sources(const BlockName& block_name) const;

    /**
     * @brief Returns the source parameters of a given parameter.
     *
     * For plain @ref Parameter objects, this is usually an empty map.
     * For @ref DependentParameter objects, this returns the actual source parameters.
     *
     * @param block_name Name of the block containing the parameter.
     * @param id         LHA identifier of the parameter.
     * @return Map of source ParamId to Parameter pointer.
     *
     * @throws If the block or parameter does not exist.
     */
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> get_parameter_sources(const BlockName& block_name, LhaID id) const;

    /**
     * @brief Recursively finds all leaf source parameters behind a set of parameters.
     *
     * Starting from each input @ref ParamId, this method recursively explores:
     * - parameter-level dependencies via `Parameter::get_source_parameters()`,
     * - block-level dependencies via `Block::get_source_blocks()`.
     *
     * A parameter is considered a leaf if:
     * - it cannot be resolved in this accessor, or
     * - it has no parameter/block sources.
     *
     * @param param_ids Initial set of parameters to inspect.
     * @return Set of leaf/root source parameters.
     */
    std::unordered_set<ParamId>
    get_all_source_parameters(const std::unordered_set<ParamId>& param_ids) const;
    
    /**
     * @brief Alias-aware mutable access to a block.
     *
     * @param block_name Name or alias of the block.
     * @return Reference to the stored block shared pointer.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    std::shared_ptr<Block>& at(const BlockName& block_name);

    /**
     * @brief Alias-aware const access to a block.
     *
     * @param block_name Name or alias of the block.
     * @return Const reference to the stored block shared pointer.
     *
     * @throws std::invalid_argument if the block cannot be resolved.
     */
    const std::shared_ptr<Block>& at(const BlockName& block_name) const;

    /**
     * @brief Merges two block accessors with no conflict allowed.
     *
     * All blocks from @p lhs and @p rhs are deep-copied into a new accessor
     * using the `Block(std::shared_ptr<Block>)` constructor.
     *
     * If both accessors contain a block with the same resolved name, this is
     * considered an error.
     *
     * @param lhs Left-hand accessor.
     * @param rhs Right-hand accessor.
     * @return New accessor containing all blocks from both operands.
     *
     * @throws Via logging/error handling if overlapping block names are found.
     */
    friend std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Merges two accessors with priority given to @p lhs.
     *
     * Semantics:
     * - all rhs blocks are copied first,
     * - lhs blocks are then merged on top,
     * - if a block exists only in lhs, it is copied entirely,
     * - if a parameter exists in both lhs and rhs:
     *   - lhs central value is kept,
     *   - if lhs uncertainties are non-zero, lhs is kept unchanged,
     *   - if lhs has zero stat+syst uncertainties, rhs uncertainties are copied into
     *     a copy of lhs.
     *
     * @param lhs Priority accessor.
     * @param rhs Secondary accessor.
     * @return New merged accessor.
     */
    friend std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Extracts a subset of blocks into a new accessor.
     *
     * This operation performs a shallow copy of block shared pointers:
     * the returned accessor points to the same underlying @ref Block objects.
     *
     * @param block_names Set of block names to extract.
     * @return New accessor containing the requested subset.
     *
     * @throws If one requested block does not exist.
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
     * @brief Inserts or replaces a block, updating alias metadata.
     *
     * The block is stored under a canonical internal key chosen from its aliases
     * (or merged into an existing key if one alias is already known).
     *
     * After insertion, `Block::bind_self()` is called on the stored block so it
     * can later retrieve a valid weak self reference.
     *
     * @param name Logical block name (possibly with several aliases).
     * @param blk  Block instance to store.
     */
    void emplace(const BlockName& name, std::shared_ptr<Block> blk);
    
    /**
     * @brief Detaches a dependent block from its upstream dependencies.
     *
     * This method resolves the block, checks that it is actually a
     * @ref DependentBlock, and then calls `DependentBlock::detach()`.
     *
     * Detaching keeps the current cached contents of the block but removes the
     * active dependency links to source blocks until reattachment.
     *
     * @param block_name Name of the block to detach.
     * If the resolved block is not a DependentBlock, do nothing.
     */
    void detach_block(const BlockName& block_name);

    /**
     * @brief Reattaches a previously detached dependent block to its saved sources.
     *
     * This method resolves the block, checks that it is a @ref DependentBlock,
     * and then calls `DependentBlock::reattach()`.
     *
     * @param block_name Name of the block to reattach.
     * If the resolved block is not a DependentBlock, do nothing.
     */

    void reattach_block(const BlockName& block_name);

    /**
     * @brief Detaches a dependent parameter from its upstream dependencies.
     *
     * The target parameter is retrieved from the given block, dynamically cast to
     * @ref DependentParameter, and then `DependentParameter::detach()` is called.
     *
     * @param block_name Name of the block containing the parameter.
     * @param id         LHA identifier of the parameter to detach.
     * If the resolved parameter is not a DependentParameter, do nothing.
     */
    void detach_parameter(const BlockName& block_name, LhaID id);

    /**
     * @brief Reattaches a previously detached dependent parameter.
     *
     * The target parameter is retrieved from the given block, dynamically cast to
     * @ref DependentParameter, and then `DependentParameter::reattach()` is called.
     *
     * @param block_name Name of the block containing the parameter.
     * @param id         LHA identifier of the parameter to reattach.
     * If the resolved parameter is not a DependentParameter, do nothing.
     */
    void reattach_parameter(const BlockName& block_name, LhaID id);

    /**
     * @brief Erases a block and removes all associated aliases from metadata.
     *
     * If the name cannot be resolved, this is a no-op.
     *
     * @param name Block name or alias.
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
    
    /// Alias-aware scale existence check on a raw name/alias.
    bool has_scale(std::string_view block) const;

    /// Alias-aware scale getter on a raw name/alias.
    double get_scale(std::string_view block) const;
    
    /// Alias-aware existence check on a raw name/alias.
    bool contains(std::string_view name) const;

    /// Convenience insertion overload from raw string view.
    void emplace(std::string_view name, std::shared_ptr<Block> blk) {
        emplace(BlockName(std::string(name)), std::move(blk));
    }

    /// Alias-aware mutable lookup by raw name/alias.
    std::shared_ptr<Block>& at(std::string_view name);

    /// Alias-aware const lookup by raw name/alias.
    const std::shared_ptr<Block>& at(std::string_view name) const;
    
    /// Alias-aware parameter existence check.
    bool has_param(std::string_view block, const LhaID& id) const;

    /// Alias-aware value getter.
    scalar_t getValue(std::string_view block, const LhaID& id) const;

    /// Alias-aware value setter.
    void setValue(std::string_view block, const LhaID& id, scalar_t value);

    /// Alias-aware parameter removal.
    void remove_item(std::string_view block, const LhaID& id);

    /// Alias-aware block erasure by name/key.
    void erase_block(std::string_view alias_or_key);
    
    /**
     * @brief Returns the normalized alias set of a BlockName.
     *
     * Normalization is uppercase-based through @ref normalize().
     */
    static std::unordered_set<std::string> norm_aliases(const BlockName& n);

    /**
     * @brief Returns the normalized alias set of a BlockName.
     *
     * Normalization is uppercase-based through @ref normalize().
     */
    static std::string choose_key(const std::unordered_set<std::string>& aliases_norm);

    /**
     * @brief Normalizes one block alias.
     *
     * Current normalization policy converts characters to uppercase.
     *
     * @param s Raw alias.
     * @return Normalized alias.
     */
    static std::string normalize(std::string_view s);
    
    /**
     * @brief Returns the canonical key corresponding to an alias, if known.
     *
     * @param alias Alias/name candidate.
     * @return Canonical storage key, or empty string if unresolved.
     */
    std::string key_for(std::string_view alias) const;

    /**
     * @brief Returns the canonical key corresponding to a BlockName, if known.
     *
     * Tries all aliases carried by the BlockName.
     *
     * @param name Block name object.
     * @return Canonical storage key, or empty string if unresolved.
     */
    std::string key_for(const BlockName& name) const;

    /**
     * @brief Returns the canonical key corresponding to a std::string alias.
     *
     * @param name Alias/name candidate.
     * @return Canonical storage key, or empty string if unresolved.
     */
    std::string key_for(const std::string& name) const;
    
    /**
     * @brief Resolves a raw alias/name into a canonical key.
     */
    std::string resolve_key(std::string_view name) const;

    /**
     * @brief Resolves a raw alias/name into a canonical key.
     */
    std::string resolve_key(const BlockName& name) const;

    /**
     * @brief Merges a new BlockName into the metadata of an existing canonical key.
     *
     * All aliases from @p name are added to the metadata entry stored under @p key,
     * and alias-to-key mappings are updated accordingly.
     *
     * @param key  Canonical key already present in storage.
     * @param name New block name / aliases to merge.
     */
    void merge_name_into_key(const std::string& key, const BlockName& name);
    
    /// Maps normalized aliases to canonical storage keys.
    std::unordered_map<std::string, std::string> alias_to_key_;

    /// Maps canonical storage keys back to full logical BlockName objects.
    std::unordered_map<std::string, BlockName> key_to_name_;
};

#endif