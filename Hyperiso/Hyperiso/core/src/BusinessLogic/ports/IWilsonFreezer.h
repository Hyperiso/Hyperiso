#ifndef IWILSONFREEZER_H
#define IWILSONFREEZER_H

#include "Freezer.h"

/**
 * @file IWilsonFreezer.h
 * @brief Interface for freezing/unfreezing Wilson blocks in the dependency system.
 *
 * Freezing a block typically prevents recomputation (or locks parameter values)
 * while scanning nuisance parameters, doing fits, or performing repeated observable
 * evaluations at fixed Wilson coefficients.
 *
 * Implementations decide how to map a user-level identifier (group, block name, etc.)
 * to the actual set of internal parameter blocks to freeze/unfreeze.
 *
 * @tparam T Identifier type used by the implementation (e.g. @ref WGroupId).
 *
 * @see WilsonFreezer
 * @see Freezer
 */
template<typename T>
class IWilsonFreezer {
public:
    /**
     * @brief Freeze all parameter blocks associated to @p id.
     */
    virtual void freeze(T) = 0;

    /**
     * @brief Unfreeze all parameter blocks associated to @p id.
     */
    virtual void unfreeze(T) = 0;
};

#endif // IWILSONFREEZER_H
