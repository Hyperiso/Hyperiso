#ifndef PARAMETERROUTER_H
#define PARAMETERROUTER_H

#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "Include.h"
#include "Utils.h"

/**
 * @file ParameterRouter.h
 * @brief Defines structures and classes for routing and categorizing parameter blocks.
 *
 * This file contains:
 * - ParameterBlockRepartition: a structure that maps blocks to their parameter types.
 * - ParametersAccessRights: a structure that defines access rights to parameters for various models.
 * - ParamRouter: a class that provides utility functions to determine the type of a parameter block.
 */

/**
 * @struct ParameterBlockRepartition
 * @brief Contains the repartition of parameter blocks according to their ParameterType.
 *
 * This is essentially a static "routing table" telling the framework which
 * LHA / FLHA blocks conceptually belong to each high-level parameter family
 * (SM, BSM, flavor, Wilson, decay, observable, passthrough).
 */
struct ParameterBlockRepartition {
    /**
     * @brief Static map associating a ParameterType with a set of block names.
     *
     * Interpretation:
     *   - ParameterType::SM      → blocks that are purely SM or baseline inputs.
     *   - ParameterType::BSM     → blocks containing BSM spectrum/mixing/etc.
     *   - ParameterType::FLAVOR  → FLHA flavor inputs (masses, lifetimes, etc.).
     *   - ParameterType::WILSON  → Wilson coefficients and scale tags.
     *   - ParameterType::DECAY   → decay-related internal blocks.
     *   - ParameterType::OBSERVABLE → observable blocks (FOBS, FDIPOLE, ...).
     *   - ParameterType::PASSTHROUGH → blocks simply forwarded / copied.
     *
     * Note: some blocks such as "MASS" and "GAUGE" can be used by multiple
     * parameter types (e.g. SM vs BSM). In those cases, finer logic is
     * implemented in ParamRouter::GetType(block, id).
     */
    static inline const std::map<ParameterType, std::unordered_set<BlockName>> BLOCKS {
        {ParameterType::SM, {"SMINPUTS", "MASS", "VCKMIN", "UPMNSIN", "UPMNS", "IMUPMNS", "VCKM", "GAUGE"}}, //TODO VCKM, GAUGE
        {ParameterType::BSM, {"GAUGE", "MASS", "HMIX", "ALPHA", "MSOFT", "NMIX", "UMIX", "VMIX", "NMAMIX", "NMHMIX", "STOPMIX", "SBOTMIX", "STAUMIX", "AU", "AD", "AE", "YU", "YD", "YE", "MINPAR"}},
        {ParameterType::FLAVOR, {"FMASS", "FLIFE", "FCONST", "FCONSTRATIO", "FBAG", "FPARAM"}},
        {ParameterType::WILSON, {"FWCOEF", "IMFWCOEF", "EW_SCALE", "B_SCALE", "D_SCALE", "K_SCALE"}},
        {ParameterType::DECAY, {"B_Ks", "B_ll", "B_Xs", "B_Dlnu", "B_Dslnu", "B_Xsll", "B_Ksll", "M0_Mix", "B_phi", "B_K", "K_ll", "K_pi", "K_lnu", "Lb_L"}},
        {ParameterType::OBSERVABLE, {"FOBS", "FOBSERR", "FOBSSM", "FOBSSMERR", "FDIPOLE"}},
        {ParameterType::PASSTHROUGH, {"MODSEL", "SPINFO", "FMODSEL", "FCINFO", "EXTPAR"}}
    };

    /**
     * @brief Filters custom (non-standard) blocks from a source list.
     *
     * A block is considered "custom" if its lower-cased name is not found in
     * any of the standard BLOCKS sets, with two special exceptions:
     *   - "mass"
     *   - "gauge"
     * which are always considered custom if present in @p source.
     *
     * @param source A list of input block names (as found in a file or DB).
     * @return A list of block names considered as custom blocks.
     */
    static std::vector<BlockName> filter_custom_blocks(const std::vector<BlockName>& source);
};

/**
 * @struct ParametersAccessRights
 * @brief Contains per-model access rights to parameters in shared blocks.
 *
 * Some blocks (notably MASS and GAUGE) are used by both SM and BSM sectors.
 * This structure encodes which (block, PDG-code) entries "belong" to which
 * model, so ParamRouter can resolve the type of an individual parameter
 * even when the block is ambiguous.
 */
struct ParametersAccessRights {
    /**
     * @brief Access rights for Standard Model (SM) parameters.
     *
     * For each block, we list the set of PDG codes (or indices) that are
     * considered SM-owned. Other codes in the same block may be considered
     * BSM (or something else) by ParamRouter.
     */
    static inline const std::map<BlockName, std::unordered_set<long>> SM_RIGHTS {
        {"MASS", {1, 2, 3, 4, 11, 12, 13, 14, 15, 16, 21, 22, 24, 25}}, 
        {"GAUGE", {1, 2, 3}},
    };

    /**
     * @brief Access rights for Two-Higgs-Doublet Model (THDM) parameters.
     *
     * THDM is a particular BSM extension; this map encodes which entries
     * are relevant for that model. It does not participate directly in
     * ParamRouter::GetType, but can be used by higher-level logic if needed.
     */
    static inline const std::map<BlockName, std::unordered_set<long>> THDM_RIGHTS {
        {"MASS", {25, 35, 36, 37}},
        {"ALPHA", {0}},
        {"HMIX", {}},
        {"GAUGE", {}},
    };

    /**
     * @brief Access rights for Supersymmetric (SUSY) model parameters.
     *
     * Similar to THDM_RIGHTS but for a SUSY-like model. Again, this is
     * primarily a data structure for model-level logic outside ParamRouter.
     */
    static inline const std::map<BlockName, std::unordered_set<long>> SUSY_RIGHTS {
        {"MASS", {25, 35, 36, 37, 
                  1000001, 1000002, 1000003, 1000004, 1000005, 1000006, 1000011, 1000012, 1000013, 1000014, 1000015, 1000016, 
                  2000001, 2000002, 2000003, 2000004, 2000005, 2000006, 2000011, 2000013, 2000015, 
                  1000021, 1000022, 1000023, 1000024, 1000025, 1000035, 1000037, 1000039}}, 
        {"GAUGE", {1, 2, 3}},
    };
};

/**
 * @class ParamRouter
 * @brief Provides utilities for identifying the type of parameter blocks/entries.
 *
 * The router answers questions of the form:
 *   - "Given this LHA block and ID, is this parameter SM or BSM or ...?"
 *   - "Which ParameterTypes own this block?"
 *   - "Which blocks are associated with a given ParameterType?"
 *
 * It relies on:
 *   - ParameterBlockRepartition::BLOCKS
 *   - ParametersAccessRights::SM_RIGHTS (for resolving ambiguous blocks).
 */
class ParamRouter {
public:
    /**
     * @brief Determines the ParameterType of a given (block, LhaID) pair.
     *
     * Resolution logic:
     *   1. If the block appears in ParametersAccessRights::SM_RIGHTS,
     *      then:
     *         - if the ID is in SM_RIGHTS[block] → ParameterType::SM
     *         - otherwise                       → ParameterType::BSM
     *   2. Otherwise, the function scans ParameterBlockRepartition::BLOCKS
     *      and returns the first ParameterType whose block set contains @p block.
     *   3. If no match is found, it logs an error and does not return.
     *
     * @param block Name of the block.
     * @param id    LHA ID of the parameter (used for conflict resolution).
     * @return The corresponding ParameterType.
     *
     * @note On error, this function calls LOG_ERROR, which typically terminates
     *       the program or throws, depending on your Logger implementation.
     */
    static ParameterType GetType(BlockName block, LhaID id);

    /**
     * @brief Retrieves all ParameterTypes associated with a block.
     *
     * This can be useful when a block is shared by multiple categories
     * (e.g. "MASS" for SM and BSM). Unlike the (block,id) version, this
     * does not attempt to resolve conflicts; it simply reports all types
     * that list @p block in ParameterBlockRepartition::BLOCKS.
     *
     * @param block Name of the block.
     * @return A list of ParameterTypes corresponding to the block.
     */
    static std::vector<ParameterType> GetType(BlockName block);

    /**
     * @brief Retrieves all blocks owned by a given ParameterType.
     *
     * Thin wrapper around ParameterBlockRepartition::BLOCKS.at(ptype).
     *
     * @param ptype The ParameterType of interest.
     * @return A set of block names belonging to the specified type.
     *
     * @throws std::out_of_range if @p ptype is not present in BLOCKS.
     */
    static std::unordered_set<BlockName> GetOwnedBlocks(ParameterType ptype);
   
};

#endif // PARAMETERROUTER_H
