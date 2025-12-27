#ifndef USE_MARTY_H
#define USE_MARTY_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

/**
 * @file UseMarty.h
 * @brief Core API indicating whether the MARTY backend is active.
 *
 * This header defines @ref UseMarty, a boolean-valued @ref ICoreAPI
 * that reports whether the current model is @ref Model::MARTY.
 */

/**
 * @class UseMarty
 * @ingroup CoreAPIModule
 * @brief Core API returning whether MARTY is used as backend.
 *
 * This API is typically used to conditionally enable MARTY-specific
 * logic without directly querying configuration internals.
 */
class UseMarty : public ICoreAPI<bool> {
public:
    /**
     * @brief Indicates whether the current model is MARTY.
     *
     * @return true if the active model is Model::MARTY, false otherwise.
     */
    inline bool get() override {return HyperisoMaster().get_model() == Model::MARTY ? 1 : 0;}
};

#endif // USE_MARTY_H