/**
 * @file BlockAccessor.h
 * @brief Provides an accessor for multiple parameter blocks.
 * 
 * This file defines the BlockAccessor class, which manages multiple blocks
 * and provides a unified interface for accessing and modifying parameters.
 */
#if !defined(BLOCK_ACCESSOR_H)
#define BLOCK_ACCESSOR_H

#include <map>
#include <string>
#include <memory>
#include <stdexcept>
#include "Block.h"

/**
 * @class BlockAccessor
 * @brief A class that manages multiple parameter blocks, for the Parameters class.
 * 
 * This class allows adding, retrieving, and modifying values in multiple parameter blocks.
 */
class BlockAccessor : public Block {
public:
    /**
     * @brief Adds a new block to the accessor.
     * @param name The name of the block.
     * @param block A shared pointer to the block.
     */
    void addBlock(const std::string& name, std::shared_ptr<Block> block);
    
    /**
     * @brief Checks if a block exists with a given parameter.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @return True if the block exists, false otherwise.
     */
    bool exist(const std::string blockName, LhaID pdgCode) const;

    /**
     * @brief Checks if a block exists.
     * @param blockName The name of the block.
     * @return True if the block exists, false otherwise.
     */
    bool has_block(const std::string blockName) const;

    /**
     * @brief Sets the value of a parameter in a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to set.
     * @param force If true, forces the update.
     */
    void setValue(const std::string& blockName, LhaID pdgCode, double value, bool force = false);

    /**
     * @brief Sets the mode of a parameter in a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param mode The mode to set.
     */
    void setMode(const std::string& blockName, LhaID pdgCode, ParameterMode mode);

    /**
     * @brief Retrieves the value of a parameter from a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @return The parameter value.
     */
    double getValue(const std::string& blockName, LhaID pdgCode) const;

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
    std::vector<std::string> get_block_names();

    /**
     * @brief Retrieves a block from the stored blocks.
     * @return A shared_ptr to the given block.
     */
    std::shared_ptr<Block> get_block(const std::string& block_name);

    /**
     * @brief Override to prevent direct access to all values.
     * @throws std::logic_error Always throws an error.
     */
    std::map<LhaID, double> getAllValues() override{
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct access to parameter values.
     * @throws std::logic_error Always throws an error.
     */
    double getValue(LhaID pdgCode) const override {
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct setting of parameter values.
     * @throws std::logic_error Always throws an error.
     */
    void setValue(LhaID pdgCode, double value, bool force = false) override {
        throw std::logic_error("Use setValue with block name for BlockAccessor");
    }

     /**
     * @brief Override to prevent direct setting of parameter modes.
     * @throws std::logic_error Always throws an error.
     */
    void setDeviation(LhaID id, double std_stat, double std_syst, bool force = false) {
        throw std::logic_error("Use setDeviation with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct setting of parameter modes.
     * @throws std::logic_error Always throws an error.
     */
    void setMode(LhaID pdgCode, ParameterMode mode) override {
        throw std::logic_error("Use setMode with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct access to parameters.
     * @throws std::logic_error Always throws an error.
     */
    std::vector<LhaID> getAllIDs() {
        throw std::logic_error("Use getAllIDs with block name for BlockAccessor");
    }

    /**
     * @brief Retrieves all parameter ids.
     * @return A vector of LhaIDs of stored parameters.
     */
    bool hasID(LhaID id) {
        throw std::logic_error("Use exists with block name for BlockAccessor");
    }


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

    friend std::ostream& operator<<(std::ostream&, std::shared_ptr<BlockAccessor>);

private:
    /// A map of block names to their corresponding shared pointers.
    std::map<std::string, std::shared_ptr<Block>> blocks;
};

#endif