#ifndef MARTYMODELPATHAPI_H
#define MARTYMODELPATHAPI_H

#include "MartyAdapter.h"
#include "Config.h"
#include "ICoreAPI.h"

/**
 * @file MartyModelPathAPI.h
 * @brief Core API exposing the filesystem path to the MARTY model file.
 *
 * This header defines @ref MartyModelPathAPI, a concrete implementation
 * of @ref ICoreAPI returning the path to the MARTY model definition file
 * (typically a `.h` file).
 *
 * The path is resolved through @ref MartyAdapter and ultimately comes
 * from the active @ref Config used to initialize Hyperiso.
 *
 * @see ICoreAPI
 * @see MartyAdapter
 * @see Config
 */

/**
 * @class MartyModelPathAPI
 * @ingroup CoreAPIModule
 * @brief Core API returning the path to the MARTY model file.
 *
 * This API allows clients to query where the MARTY model definition
 * file is located on disk, without directly accessing configuration
 * internals or singletons.
 */
class MartyModelPathAPI : public ICoreAPI<fs::path> {
public:
    /**
     * @brief Returns the filesystem path to the MARTY model file.
     *
     * @return Path to the MARTY model file.
     */
    inline fs::path get() override { return MartyAdapter().get_path(MartyPath::MODEL_FILE); }
}; 

#endif // MARTYMODELNAMEAPI_H
