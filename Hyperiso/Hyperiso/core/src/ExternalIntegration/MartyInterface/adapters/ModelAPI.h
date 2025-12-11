#ifndef MODEL_API_H
#define MODEL_API_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "ICoreAPI.h"

/**
 * @file ModelAPI.h
 * @brief Provides a thin ICoreAPI wrapper around HyperisoMaster::get_model().
 *
 * This header defines ::ModelAPI, an ::ICoreAPI<Model> implementation
 * that queries the currently active physics model from ::HyperisoMaster.
 */

/**
 * @class ModelAPI
 * @ingroup MonitoringSystemModule
 * @brief Core API exposing the current physics model.
 *
 * ModelAPI is a convenience class that queries ::HyperisoMaster to
 * obtain the ::Model currently used by the framework. It is typically
 * injected into components that need to know whether the working model
 * is SM, THDM, SUSY, etc., without directly depending on HyperisoMaster.
 */
class ModelAPI : public ICoreAPI<Model> {
public:
    /**
     * @brief Returns the current ::Model value.
     *
     * Internally this calls ::HyperisoMaster::get_model().
     */
    inline Model get() override {return HyperisoMaster().get_model();}
};

#endif  // MODEL_API_H