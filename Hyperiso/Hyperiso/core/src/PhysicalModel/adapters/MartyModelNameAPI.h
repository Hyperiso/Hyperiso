#ifndef MARTYMODELNAMEAPI_H
#define MARTYMODELNAMEAPI_H

#include "MartyAdapter.h"
#include "Config.h"
#include "ICoreAPI.h"

/**
 * @file MartyModelNameAPI.h
 * @brief Core API exposing the active MARTY model name.
 *
 * This header defines @ref MartyModelNameAPI, a small adapter implementing
 * @ref ICoreAPI to expose the name of the currently active MARTY model
 * as configured in the Hyperiso framework.
 *
 * The returned value corresponds to:
 * - the MARTY class name (e.g. "THDM", "MSSM", etc.),
 * - typically provided through the HyperisoConfig object and accessed via
 *   higher-level adapters.
 *
 * @see ICoreAPI
 * @see MartyAdapter
 * @see HyperisoConfig
 */

/**
 * @class MartyModelNameAPI
 * @ingroup CoreAPIModule
 * @brief Core API returning the MARTY model class name.
 *
 * This API provides read-only access to the MARTY model name currently
 * configured in Hyperiso.
 *
 * Internally, it delegates to @ref MartyAdapter, which encapsulates
 * access to MARTY-related configuration details.
 *
 * Typical usage:
 * @code
 *   std::shared_ptr<ICoreAPI<std::string>> api =
 *       std::make_shared<MartyModelNameAPI>();
 *
 *   std::string modelName = api->get();
 * @endcode
 */
class MartyModelNameAPI : public ICoreAPI<std::string> {
public:
    /**
     * @brief Returns the name of the active MARTY model.
     *
     * @return MARTY model class name as a string.
     */
    inline std::string get() override { return MartyAdapter().get_marty_model_name(); }
}; 

#endif // MARTYMODELNAMEAPI_H
