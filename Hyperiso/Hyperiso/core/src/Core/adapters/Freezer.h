#ifndef FREEZER_H
#define FREEZER_H

#include "IFreezer.h"
#include "Include.h"
#include "Parameters.h"

/**
 * @class Freezer
 * @ingroup FreezeControlModule
 * @brief Concrete implementation for freezing and unfreezing blocks and parameters.
 *
 * Provides static methods to freeze or unfreeze parameters and blocks by delegating to Parameters system.
 */
class Freezer : public IFreezer<ParameterType, std::string, ParamId> {
public:

    /**
     * @brief Freezes a block given its type and name.
     * @param p_type The parameter type (e.g., SM, BSM).
     * @param block_name The name of the block to freeze.
     */
    static void freeze(const ParameterType&, const std::string&);

    /**
     * @brief Freezes a single parameter.
     * @param pid The parameter ID to freeze.
     */
    static void freeze(const ParamId&);

    /**
     * @brief Unfreezes a block given its type and name.
     * @param p_type The parameter type.
     * @param block_name The name of the block to unfreeze.
     */
    static void unfreeze(const ParameterType&, const std::string&);

    /**
     * @brief Unfreezes a single parameter.
     * @param pid The parameter ID to unfreeze.
     */
    static void unfreeze(const ParamId&);
};

#endif // FREEZER_H
