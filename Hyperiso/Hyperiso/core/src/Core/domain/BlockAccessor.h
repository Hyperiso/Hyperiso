#ifndef BLOCK_ACCESSOR_H
#define BLOCK_ACCESSOR_H

#include <functional>
#include <unordered_map>

#include "Include.h"
#include "Block.h"

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
 */
class BlockAccessor : public std::unordered_map<BlockName, std::shared_ptr<Block>> {
public:
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
     */
    std::unordered_set<BlockName> get_block_names();

    /**
     * @brief Removes a parameter from a specified block.
     *
     * Calls Block::remove(id) on the underlying block. If the block does not
     * exist, a warning is logged and no action is taken.
     *
     * @param block_name The name of the block.
     * @param id         The LhaID of the parameter to remove.
     */
    void remove_item(const BlockName& block_name, LhaID id);

    /**
     * @brief Checks if a given block exists in the accessor.
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
     * @brief Accessor to a block by name (non-const).
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
     */
    friend std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Extracts a subset of blocks from the accessor.
     *
     * Only blocks listed in @p block_names are inserted in the result.
     * If any requested block does not exist, an error is logged.
     *
     * @param block_names A set of block names to retrieve.
     * @return A shared pointer to a new BlockAccessor containing the requested blocks.
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
};

#endif