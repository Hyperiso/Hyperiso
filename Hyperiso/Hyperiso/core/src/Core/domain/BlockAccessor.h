/**
 * @file BlockAccessor.h
 * @brief Provides an accessor for multiple parameter blocks.
 * 
 * This file defines the BlockAccessor class, which manages multiple blocks
 * and provides a unified interface for accessing and modifying parameters.
 */
#if !defined(BLOCK_ACCESSOR_H)
#define BLOCK_ACCESSOR_H

#include <functional>
#include <unordered_map>

#include "Include.h"
#include "Block.h"

/**
 * @class BlockAccessor
 * @brief A class that manages multiple parameter blocks, for the Parameters class.
 * 
 * This class allows adding, retrieving, and modifying values in multiple parameter blocks.
 */
class BlockAccessor : public std::unordered_map<BlockName, std::shared_ptr<Block>> {
public:
    /**
     * @brief Checks if a block exists with a given parameter.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @return True if the block exists, false otherwise.
     */
    bool has_param(const BlockName& blockName, LhaID pdgCode) const;

    /**
     * @brief Sets the value of a parameter in a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to set.
     */
    void setValue(const BlockName& blockName, LhaID pdgCode, scalar_t value);

    // /**
    //  * @brief Sets the mode of a parameter in a specified block.
    //  * @param blockName The name of the block.
    //  * @param pdgCode The PDG code of the parameter.
    //  * @param mode The mode to set.
    //  */
    // void setMode(const std::string& blockName, LhaID pdgCode, ParameterMode mode);

    /**
     * @brief Retrieves the value of a parameter from a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @return The parameter value.
     */
    scalar_t getValue(const BlockName& blockName, LhaID pdgCode) const;

    /**
     * @brief Retrieves the parameter from a specified block.
     * @param blockName The name of the block.
     * @param id The LHA ID of the parameter.
     * @return The parameter.
     */
    std::shared_ptr<Parameter> getParameter(const BlockName& blockName, LhaID id) const;

    /**
     * @brief Retrieves the parameter from a specified block.
     * @param blockName The name of the block.
     * @param id The LHA ID of the parameter.
     * @param source The source parameter.
     */
    void setParameter(const BlockName& blockName, LhaID id, std::shared_ptr<Parameter> source);

    /**
     * @brief Retrieves all values from a specified block.
     * @param blockName The name of the block.
     * @return A map of PDG codes to parameter values.
     */
    std::map<LhaID, double> getAllValues(BlockName blockName);

    /**
     * @brief Retrieves a list of all block names.
     * @return A vector of block names.
     */
    std::unordered_set<BlockName> get_block_names();

    /**
     * @brief Removes a parameter from a specified block.
     * 
     * If the block does not exist, a warning is logged.
     * 
     * @param block_name The name of the block.
     * @param id The LHA ID of the parameter to remove.
     */
    void remove_item(const BlockName& block_name, LhaID id);

    bool contains(const BlockName& block_name) const;

    std::shared_ptr<Block>& at(const BlockName& block_name);
    const std::shared_ptr<Block>& at(const BlockName& block_name) const;

    /**
     * @brief Merges two BlockAccessor instances without resolving conflicts.
     * 
     * Throws an error if both accessors contain a block with the same name.
     *
     * @param lhs First BlockAccessor instance.
     * @param rhs Second BlockAccessor instance.
     * @return A shared pointer to a new BlockAccessor containing all unique blocks.
     * @throws Error if overlapping blocks are found.
     */
    friend std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Merges two BlockAccessor instances with priority given to lhs.
     *
     * If both accessors contain the same block, parameters from lhs overwrite those from rhs.
     *
     * @param lhs First BlockAccessor instance (priority).
     * @param rhs Second BlockAccessor instance.
     * @return A shared pointer to a new BlockAccessor containing merged blocks.
     */
    friend std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Retrieves specific blocks from the accessor.
     * 
     * Only blocks listed in `block_names` will be extracted.
     * If any block does not exist, an error is logged.
     * 
     * @param block_names A set of block names to retrieve.
     * @return A shared pointer to a new BlockAccessor containing the requested blocks.
     */
    std::shared_ptr<BlockAccessor> operator[](std::unordered_set<BlockName> block_names);
    /**
     * @brief Stream output operator for BlockAccessor.
     *
     * Prints all block names and their parameter values.
     *
     * @param os Output stream.
     * @param ba Shared pointer to the BlockAccessor to print.
     * @return The output stream.
     */
    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<BlockAccessor>);
};

#endif