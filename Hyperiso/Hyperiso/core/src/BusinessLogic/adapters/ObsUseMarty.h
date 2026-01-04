#ifndef OBS_USE_MARTY_H
#define OBS_USE_MARTY_H

#include "HyperisoMaster.h"
#include "Config.h"
#include "IObsCoreAPI.h"

/**
 * @file ObsUseMarty.h
 * @brief Observable-layer API exposing whether MARTY backend is active.
 *
 * This small adapter lets the observable/business layer query the global runtime
 * configuration to know if the current model backend is MARTY.
 *
 * It is intentionally read-only and returns a simple boolean.
 *
 * @see HyperisoMaster
 * @see Model
 * @see IObsCoreAPI
 */
class ObsUseMarty : public IObsCoreAPI<bool> {
public:
    /**
     * @brief Returns true if the selected model is @ref Model::MARTY.
     */
    inline bool get() override {return HyperisoMaster().get_model() == Model::MARTY ? 1 : 0;}
};

#endif // OBS_USE_MARTY_H