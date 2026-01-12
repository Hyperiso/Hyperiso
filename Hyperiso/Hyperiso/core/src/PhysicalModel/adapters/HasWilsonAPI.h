#ifndef HasWilsonAPI_H
#define HasWilsonAPI_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

/**
 * @file HasWilsonAPI.h
 * @brief Core API exposing whether Wilson coefficients are provided as user inputs.
 *
 * This header defines @ref HasWilsonAPI, a concrete @ref ICoreAPI<bool> implementation.
 * It queries Hyperiso's active configuration flag @ref ExternalFlag::HAS_WILSON_INPUT
 * through @ref HyperisoMaster.
 *
 * This is typically used by higher-level components that need to adjust their behavior
 * when Wilson coefficients are treated as external inputs (e.g. preventing modifications
 * to certain matching scales, or changing initialization flows).
 *
 * @see HyperisoMaster
 * @see HyperisoConfig
 * @see ExternalFlag
 */

/**
 * @class HasWilsonAPI
 * @ingroup MonitoringSystemModule
 * @brief Returns whether Hyperiso currently uses Wilson coefficients as external inputs.
 *
 * This class is a lightweight adapter over @ref HyperisoMaster:
 * @code
 *   HasWilsonAPI hasWilson;
 *   if (hasWilson.get()) { ... }
 * @endcode
 *
 * Internally, it delegates to:
 * @code
 *   HyperisoMaster().check_flag(ExternalFlag::HAS_WILSON_INPUT);
 * @endcode
 *
 * @note This reads the flag from the active HyperisoConfig (stored in the framework cache).
 */
class HasWilsonAPI : public ICoreAPI<bool> {
public:
    /**
     * @brief Returns true if @ref ExternalFlag::HAS_WILSON_INPUT is active.
     *
     * @return true if Wilson coefficients are user-provided inputs, false otherwise.
     */
    inline bool get() override {return HyperisoMaster().check_flag(ExternalFlag::HAS_WILSON_INPUT);}
};

#endif