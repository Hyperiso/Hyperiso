#ifndef BLOCK_PROVIDER_H
#define BLOCK_PROVIDER_H

#include "Include.h"
#include "Parameters.h"
#include "IBlockProvider.h"

/**
 * @file BlockProvider.h
 * @brief Concrete implementation of IBlockProvider for HyperISO parameter blocks.
 *
 * This file defines BlockProvider, which uses Parameters to:
 * - check the existence of blocks for a given ParameterType,
 * - log all blocks of a given ParameterType via the logging system,
 * - log a specific block of a given ParameterType.
 */

 /**
 * @class BlockProvider
 * @brief Block provider bound to HyperISO Parameters / ParameterType.
 *
 * This class specializes IBlockProvider with:
 * - T  = ParameterType
 * - U  = const std::string&
 *
 * It delegates its operations to the corresponding Parameters instance:
 * - exists(...) queries Parameters::GetInstance(type)->get_blocks_list()
 * - log_all_blocks(...) logs the entire Parameters instance
 * - log_block(...) prints the content of a single block
 *
 * @see Parameters
 * @see ParameterType
 */
class BlockProvider : public IBlockProvider<ParameterType, const std::string&> {
public:
    /**
     * @brief Checks whether a block exists for the given ParameterType.
     *
     * @param blockname Name of the block to check.
     * @param pt ParameterType for which to check the block.
     * @return true if the block exists, false otherwise.
     */
    bool exists(const std::string& blockname, ParameterType) override;

    /**
     * @brief Logs all blocks for a given ParameterType.
     *
     * Internally calls:
     * @code
     * LOG_INFO(Parameters::GetInstance(type));
     * @endcode
     *
     * @param type ParameterType whose blocks should be logged.
     */
    void log_all_blocks(ParameterType type) override;

    /**
     * @brief Logs the content of a single block for a given ParameterType.
     *
     * Internally calls Parameters::GetInstance(type)->print_block(blockname).
     *
     * @param type ParameterType owning the block.
     * @param blockname Name of the block to log.
     */
    void log_block(ParameterType type, const std::string& blockname) override;

    /**
     * @brief Retrieve the content of a single block for a given ParameterType.
     *
     * Internally calls Parameters::GetInstance(type)->get_block_infos(blockname).
     *
     * @param type ParameterType owning the block.
     * @param blockname Name of the block to retrieve.
     * @return std::map<LhaID, scalar_t>, the block content
     */
    std::map<LhaID, scalar_t> get_block(ParameterType type, const std::string& blockname) override;
    
};

#endif