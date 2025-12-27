#ifndef MODEL_API_H
#define MODEL_API_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

/**
 * @file ModelAPI.h
 * @brief Core API exposing the currently active physics model.
 *
 * This header defines @ref ModelAPI, a simple implementation of
 * @ref ICoreAPI returning the active @ref Model enum value.
 *
 * @see HyperisoMaster
 * @see Config
 * @see Model
 */

/**
 * @class ModelAPI
 * @ingroup CoreAPIModule
 * @brief Core API returning the current Model.
 *
 * This API provides read-only access to the model selected during
 * Hyperiso initialization (e.g. SM, SUSY, THDM, MARTY).
 */
class ModelAPI : public ICoreAPI<Model> {
public:
    /**
     * @brief Returns the active physics model.
     *
     * @return Current model.
     */
    inline Model get() override {return HyperisoMaster().get_model();}
};

#endif // MODEL_API_H