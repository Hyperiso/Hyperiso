#ifndef FREEZER_H
#define FREEZER_H

#include "IFreezer.h"
#include "Include.h"
#include "Parameters.h"

/**
 * @file Freezer.h
 * @brief Concrete helper for freezing and unfreezing parameter blocks and parameters.
 *
 * This header defines the Freezer class, a concrete implementation of the
 * IFreezer interface for:
 *  - freezing / unfreezing entire blocks (by ParameterType and block name),
 *  - freezing / unfreezing individual parameters (by ParamId).
 *
 * Internally, it delegates to the Parameters system, which manages the
 * underlying Block / Parameter objects.
 */


/**
 * @class Freezer
 * @ingroup FreezeControlModule
 * @brief Concrete implementation for freezing and unfreezing blocks and parameters.
 *
 * This class provides a static API to:
 *  - freeze or unfreeze a whole block given its ParameterType and BlockName,
 *  - freeze or unfreeze a single parameter given its ParamId.
 *
 * All operations are delegated to the corresponding Parameters instance,
 * which in turn propagates the freeze/unfreeze requests to Block and
 * Parameter objects (via Block::freeze, Block::unfreeze, Parameter::freeze,
 * Parameter::unfreeze).
 *
 * @see Parameters
 * @see Block
 * @see Parameter
 */
class Freezer : public IFreezer<ParameterType, std::string, ParamId> {
public:

    /**
     * @brief Freezes a block given its type and name.
     *
     * Delegates to Parameters::GetInstance(p_type)->freeze_block(block_name).
     *
     * @param p_type The parameter type (e.g., SM, BSM).
     * @param block_name The name of the block to freeze.
     */
    static void freeze(const ParameterType&, const std::string&);

    /**
     * @brief Freezes a single parameter.
     *
     * Delegates to the appropriate Parameters instance based on pid.type
     * and calls freeze_param(pid.block, pid.code).
     *
     * @param pid The parameter ID to freeze.
     */
    static void freeze(const ParamId&);

    /**
     * @brief Unfreezes a block given its type and name.
     *
     * Delegates to Parameters::GetInstance(p_type)->unfreeze_block(block_name).
     *
     * @param p_type The parameter type.
     * @param block_name The name of the block to unfreeze.
     */
    static void unfreeze(const ParameterType&, const std::string&);

    /**
     * @brief Unfreezes a single parameter.
     *
     * Delegates to the appropriate Parameters instance based on pid.type
     * and calls unfreeze_param(pid.block, pid.code).
     *
     * @param pid The parameter ID to unfreeze.
     */
    static void unfreeze(const ParamId&);
};

#endif // FREEZER_H
