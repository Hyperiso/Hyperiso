/**
 * @file ParameterRouter.h
 * @brief Defines structures and classes for routing and categorizing parameter blocks.
 *
 * This file contains:
 * - ParameterBlockRepartition: a structure that maps blocks to their parameter types.
 * - ParametersAccessRights: a structure that defines access rights to parameters for various models.
 * - ParamRouter: a class that provides utility functions to determine the type of a parameter block.
 */
#ifndef __PARAMETERROUTER_H__
#define __PARAMETERROUTER_H__

#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "General.h"
#include "MemoryManager.h"
#include "Utils.h"

/**
 * @struct ParameterBlockRepartition
 * @brief Contains the repartition of parameter blocks according to their parameter type.
 */
struct ParameterBlockRepartition {
    /**
     * @brief Static map associating a ParameterType with a set of block names.
     */
    static inline const std::map<ParameterType, std::unordered_set<BlockName>> BLOCKS {
        {ParameterType::SM, {"SMINPUTS", "MASS", "VCKMIN", "UPMNSIN", "UPMNS", "IMUPMNS", "VCKM", "GAUGE"}}, //TODO VCKM, GAUGE
        {ParameterType::BSM, {"MASS", "HMIX", "ALPHA", "MSOFT", "NMIX", "UMIX", "VMIX", "A0MIX", "H0MIX", "STOPMIX", "SBOTMIX", "STAUMIX", "AU", "AD", "AE", "YU", "YD", "YE",  "UCOUPL", "DCOUPL", "LCOUPL"}},
        {ParameterType::FLAVOR, {"FMASS", "FLIFE", "FCONST", "FCONSTRATIO", "FBAG", "FPARAM"}},
        {ParameterType::WILSON, {"FWCOEF", "IMFWCOEF", "EW_SCALE", "B_SCALE"}},
        {ParameterType::DECAY, {"B_Ks", "B_ll", "B_Xs", "B_Dlnu", "B_Dslnu"}},
        {ParameterType::OBSERVABLE, {"FOBS", "FOBSERR", "FOBSSM", "FOBSSMERR", "FDIPOLE"}},
        {ParameterType::PASSTHROUGH, {"MODSEL", "SPINFO", "FMODSEL", "FCINFO", "MINPAR", "EXTPAR"}}
    };

    /**
     * @brief Filters custom (non-standard) blocks from a source list.
     *
     * @param source A list of input block names.
     * @return A list of block names considered as custom blocks.
     */
    static std::vector<BlockName> filter_custom_blocks(const std::vector<BlockName>& source);
};

/**
 * @struct ParametersAccessRights
 * @brief Contains access rights to parameters for different theoretical models (SM, THDM, SUSY).
 */
struct ParametersAccessRights {
    /**
     * @brief Access rights for Standard Model (SM) parameters.
     */
    static inline const std::map<BlockName, std::unordered_set<long>> SM_RIGHTS {
        {"MASS", {1, 2, 3, 4, 11, 12, 13, 14, 15, 16, 21, 22, 24, 25}}, 
        {"GAUGE", {1, 2, 3}},
    };

    /**
     * @brief Access rights for Two-Higgs-Doublet Model (THDM) parameters.
     */
    static inline const std::map<BlockName, std::unordered_set<long>> THDM_RIGHTS {
        {"MASS", {25, 35, 36, 37}},
        {"ALPHA", {0}},
        {"HMIX", {}},
        {"GAUGE", {{}}},
    };

    /**
     * @brief Access rights for Supersymmetric (SUSY) model parameters.
     */
    static inline const std::map<BlockName, std::unordered_set<long>> SUSY_RIGHTS {
        {"MASS", {25, 35, 36, 37, 
                  1000001, 1000002, 1000003, 1000004, 1000005, 1000006, 1000011, 1000012, 1000013, 1000014, 1000015, 1000016, 
                  2000001, 2000002, 2000003, 2000004, 2000005, 2000006, 2000011, 2000013, 2000015, 
                  1000021, 1000022, 1000023, 1000024, 1000025, 1000035, 1000037, 1000039}}, 
        {"GAUGE", {{}}},
    };
};

/**
 * @class ParamRouter
 * @brief Provides utilities for identifying the type of parameter blocks.
 */
class ParamRouter {
public:
    /**
     * @brief Determines the parameter type of a given block and ID.
     *
     * @param block Name of the block.
     * @param id LHA ID of the parameter.
     * @return The corresponding ParameterType.
     * @throws Logs an error if the block or ID is invalid.
     */
    static ParameterType GetType(BlockName block, LhaID id);

    /**
     * @brief Retrieves all parameter types associated with a block.
     *
     * @param block Name of the block.
     * @return A list of ParameterTypes corresponding to the block.
     */
    static std::vector<ParameterType> GetType(BlockName block);

    /**
     * @brief Retrieves all blocks owned by a given parameter type.
     *
     * @param ptype The ParameterType.
     * @return A set of block names belonging to the specified type.
     */
    static std::unordered_set<BlockName> GetOwnedBlocks(ParameterType ptype);
   
};

#endif // __PARAMETERROUTER_H__
