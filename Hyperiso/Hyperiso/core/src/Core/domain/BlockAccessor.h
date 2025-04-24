/**
 * @file BlockAccessor.h
 * @brief Provides an accessor for multiple parameter blocks.
 * 
 * This file defines the BlockAccessor class, which manages multiple blocks
 * and provides a unified interface for accessing and modifying parameters.
 */
#if !defined(BLOCK_ACCESSOR_H)
#define BLOCK_ACCESSOR_H

#include "Include.h"
#include "Block.h"
#include <functional>
#include <unordered_map>

/**
 * @class BlockAccessor
 * @brief A class that manages multiple parameter blocks, for the Parameters class.
 * 
 * This class allows adding, retrieving, and modifying values in multiple parameter blocks.
 */
class BlockAccessor : public std::unordered_map<std::string, std::shared_ptr<Block>> {
public:
    /**
     * @brief Checks if a block exists with a given parameter.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @return True if the block exists, false otherwise.
     */
    bool has_param(const std::string blockName, LhaID pdgCode) const;

    /**
     * @brief Sets the value of a parameter in a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to set.
     */
    void setValue(const std::string& blockName, LhaID pdgCode, double value);

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
    scalar_t getValue(const std::string& blockName, LhaID pdgCode) const;

    /**
     * @brief Retrieves the parameter from a specified block.
     * @param blockName The name of the block.
     * @param id The LHA ID of the parameter.
     * @return The parameter.
     */
    std::shared_ptr<Parameter> getParameter(const std::string& blockName, LhaID id) const;

    /**
     * @brief Retrieves the parameter from a specified block.
     * @param blockName The name of the block.
     * @param id The LHA ID of the parameter.
     * @param source The source parameter.
     */
    void setParameter(const std::string& blockName, LhaID id, std::shared_ptr<Parameter> source);

    /**
     * @brief Retrieves all values from a specified block.
     * @param blockName The name of the block.
     * @return A map of PDG codes to parameter values.
     */
    std::map<LhaID, double> getAllValues(std::string blockName);

    /**
     * @brief Retrieves a list of all block names.
     * @return A vector of block names.
     */
    std::unordered_set<std::string> get_block_names();

    void remove_item(const std::string& block_name, LhaID id);

    /**
     * @brief Merges two distinct BlockAccessor instances without priority. Throws an error if both have blocks in common.  
     * @param lhs First BlockAccessor instance.
     * @param rhs Second BlockAccessor instance.
     * @return A new `BlockAccessor` instance storing all the blocks from `rhs` and `lhs`. 
     */
    friend std::shared_ptr<BlockAccessor> operator+(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Merges two overlapping BlockAccessor instances with left-priority. If rhs has elements in common with lhs, they will be overwritten.
     * @param lhs First BlockAccessor instance.
     * @param rhs Second BlockAccessor instance.
     * @return A new `BlockAccessor` instance storing all the blocks from `rhs` and `lhs`.   
     */
    friend std::shared_ptr<BlockAccessor> operator>>(std::shared_ptr<BlockAccessor> lhs, std::shared_ptr<BlockAccessor> rhs);

    /**
     * @brief Retrieves given blocks from stored blocks if they exist.
     * @param block_names List of block names to retrieve.
     * @return A new `BlockAccessor` instance storing all the given blocks.   
     */
    std::shared_ptr<BlockAccessor> operator[](std::unordered_set<std::string> block_names);

    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<BlockAccessor>);
};

#endif