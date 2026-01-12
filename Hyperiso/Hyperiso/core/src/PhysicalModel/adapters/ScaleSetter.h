#ifndef SCALESETTER_H
#define SCALESETTER_H

#include "IParamSetter.h"
#include "Include.h"

/**
 * @file ScaleSetter.h
 * @brief Concrete setter for Wilson scale parameters.
 *
 * This header defines @ref ScaleSetter, an implementation of @ref IParamSetter
 * specialized for @ref ScaleType (matching scales, etc.).
 *
 * It writes into the WILSON Parameters instance using @ref ParameterSetter:
 * - ParameterType is forced to ParameterType::WILSON
 * - Block name is derived from ScaleType via ScaleTypeMapper::block(...)
 * - Entry id is fixed to 1 (convention used by scale blocks)
 *
 * @see IParamSetter
 * @see ParameterSetter
 * @see ScaleTypeMapper
 */

/**
 * @class ScaleSetter
 * @ingroup DataMutationModule
 * @brief Sets Wilson scale values through the Parameters mutation layer.
 *
 * This class provides a small “command-style” interface to set a given scale.
 *
 * Example:
 * @code
 *   ScaleSetter s(ScaleType::B_SCALE);
 *   s.set(4.2);
 *   s.switch_param(ScaleType::D_SCALE);
 *   s.set(1.3);
 * @endcode
 */
class ScaleSetter : public IParamSetter<ScaleType> {
public:
    /**
     * @brief Constructs a ScaleSetter targeting an initial scale type.
     *
     * @param scale_type Initial scale selector.
     */
    ScaleSetter(ScaleType scale_type) {this->param = scale_type;}

    /**
     * @brief Changes the targeted scale parameter.
     *
     * @param scale_type New scale selector.
     */
    void switch_param(ScaleType scale_type) override {this->param = scale_type;}

    /**
     * @brief Sets the value of the currently selected scale.
     *
     * Builds a ParamId:
     * - type  = ParameterType::WILSON
     * - block = ScaleTypeMapper::block(param)
     * - code  = 1
     *
     * Then delegates to @ref ParameterSetter::mutate.
     *
     * @param value Scale value to set.
     */
    void set(double value) override;

};

#endif // SCALESETTER_H
