/**
 * @file DependentBlockManager.h
 * @brief Provides static utilities to create, manage, and update dependent blocks and dependent parameters.
 *
 * This file defines the DependentBlockManager class, which offers functions to automate the creation of
 * DependentBlock and DependentParameter instances, mainly interacting with the Parameters singleton.
 */
#ifndef DEPENDENT_BLOCK_MANAGER_H
#define DEPENDENT_BLOCK_MANAGER_H

#include <map>
#include <unordered_map>
#include <memory>
#include <vector>
#include <iostream>
#include "DependentParameter.h"
#include "Block.h"

/**
 * @class DependentBlockManager
 * @brief Static class for handling dependent blocks and parameters across different Parameters instances.
 *
 * DependentBlockManager manages the lifecycle of blocks or parameters that are computed from others
 * and automatically updated when their sources change.
 */
class DependentBlockManager {
public:
    /**
     * @brief Adds a DependentBlock with multiple sources coming from various Parameters instances.
     *
     * Creates a new DependentBlock based on the provided sources and registers it in the destination Parameters instance.
     *
     * @param name Name of the DependentBlock to create.
     * @param source_names Map associating each source ParameterType to a list of block names.
     * @param dest Destination ParameterType where the DependentBlock will be registered.
     * @param recalculateFunc Lambda function used to recalculate the content of the DependentBlock.
     */
    static void addDependentBlock(
        std::string name,
        std::unordered_map<ParameterType, std::vector<std::string>> source_names,
        ParameterType dest,
        DepUpdateFunc recalculateFunc
    );

    /**
     * @brief Adds a DependentParameter based on multiple source parameters.
     *
     * Creates a DependentParameter whose value depends on other parameters, and stores it in the appropriate block.
     *
     * @param pid ID of the DependentParameter to create.
     * @param source_pids Set of IDs corresponding to the source parameters.
     * @param recalculateFunc Lambda function used to recalculate the value of the DependentParameter.
     */
    static void addDependentParameter(
        ParamId pid,
        std::unordered_set<ParamId> source_pids,
        DepParamUpdateFunc recalculateFunc
    );

    /**
     * @brief Removes a DependentBlock from a specified Parameters instance.
     *
     * @param name Name of the DependentBlock to remove.
     * @param src The ParameterType (model) where the block is registered.
     */
    static void removeDependentBlock(const std::string& name, ParameterType src);
    
    /**
     * @brief Manually triggers an update on a DependentBlock.
     *
     * @param name Name of the block to update.
     * @param src The ParameterType (model) where the block is located.
     */
    static void update(const std::string& name, ParameterType src);

private:
    /**
     * @brief Deleted default constructor. DependentBlockManager is a purely static class.
     */
    DependentBlockManager() = default;

    /**
     * @brief Deleted destructor. No instances of DependentBlockManager are allowed.
     */
    ~DependentBlockManager() = default;

    /**
     * @brief Deleted copy constructor.
     */
    DependentBlockManager(const DependentBlockManager&) = delete;

    /**
     * @brief Deleted copy assignment operator.
     */
    DependentBlockManager& operator=(const DependentBlockManager&) = delete;
};

#endif // DEPENDENT_BLOCK_MANAGER_H