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
 * @enum FlavorParamType
 * @brief Enumeration for different flavor parameter types.
 */
enum class FlavorParamType {
    LIFETIME,               ///< Lifetime parameter
    DECAY_CONSTANT,         ///< Decay constant parameter
    DECAY_CONSTANT_RATIO    ///< Ratio of decay constants
};

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
    bool exist(const std::string blockName, int pdgCode) const;

    /**
     * @brief Sets the value of a parameter in a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to set.
     * @param force If true, forces the update.
     */
    void setValue(const std::string& blockName, int pdgCode, double value, bool force = false);

    /**
     * @brief Sets the mode of a parameter in a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @param mode The mode to set.
     */
    void setMode(const std::string& blockName, int pdgCode, ParameterMode mode);

    /**
     * @brief Retrieves the value of a parameter from a specified block.
     * @param blockName The name of the block.
     * @param pdgCode The PDG code of the parameter.
     * @return The parameter value.
     */
    double getValue(const std::string& blockName, int pdgCode) const;

    /**
     * @brief Retrieves all values from a specified block.
     * @param blockName The name of the block.
     * @return A map of PDG codes to parameter values.
     */
    std::map<int, double> getAllValues(std::string blockName);

    /**
     * @brief Retrieves a list of all block names.
     * @return A vector of block names.
     */
    std::vector<std::string> get_blocks();

    /**
     * @brief Override to prevent direct access to all values.
     * @throws std::logic_error Always throws an error.
     */
    std::map<int, double> getAllValues() override{
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct access to parameter values.
     * @throws std::logic_error Always throws an error.
     */
    double getValue(int pdgCode) const override {
        throw std::logic_error("Use getValue with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct setting of parameter values.
     * @throws std::logic_error Always throws an error.
     */
    void setValue(int pdgCode, double value, bool force = false) override {
        throw std::logic_error("Use setValue with block name for BlockAccessor");
    }

    /**
     * @brief Override to prevent direct setting of parameter modes.
     * @throws std::logic_error Always throws an error.
     */
    void setMode(int pdgCode, ParameterMode mode) override {
        throw std::logic_error("Use setMode with block name for BlockAccessor");
    }

private:
    /// A map of block names to their corresponding shared pointers.
    std::map<std::string, std::shared_ptr<Block>> blocks;
};

#endif