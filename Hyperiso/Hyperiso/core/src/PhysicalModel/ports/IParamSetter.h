#ifndef IPARAMSETTER_H
#define IPARAMSETTER_H

#include "ParameterSetter.h"

/**
 * @file IParamSetter.h
 * @brief Interface for setting a specific “current parameter target” (typically a scale).
 *
 * This header defines @ref IParamSetter, an abstraction used when a component
 * conceptually manipulates a *selected parameter* (e.g. a scale parameter),
 * and needs to:
 *  - switch which parameter is targeted,
 *  - set its value.
 *
 * Typical concrete implementation:
 * - @ref ScaleSetter, which writes into ParameterType::WILSON at a block determined
 *   by the chosen ScaleType.
 *
 * @tparam T Type identifying which parameter is targeted (e.g. ScaleType).
 */

/**
 * @class IParamSetter
 * @ingroup DataMutationModule
 * @brief Abstract interface for setting a chosen parameter.
 *
 * This interface stores the currently selected target identifier in @ref param.
 *
 * @tparam T Target selector type.
 */
template<typename T>
class IParamSetter {
public:
    virtual ~IParamSetter() = default;

    /**
     * @brief Sets the numeric value of the currently selected target parameter.
     *
     * @param value Value to set.
     */
    virtual void set(double value) = 0;

    /**
     * @brief Switches which target parameter is controlled by this setter.
     *
     * @param param New target selector.
     */
    virtual void switch_param(T scale_type) = 0;

protected:
    /// Current target selector.
    T param;
};

#endif // IPARAMSETTER_H
