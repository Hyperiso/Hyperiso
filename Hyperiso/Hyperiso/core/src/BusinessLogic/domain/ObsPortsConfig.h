#ifndef OBS_PORTS_CONFIG_H
#define OBS_PORTS_CONFIG_H

#include <memory>
#include <utility>

#include "IObsParameterProxy.h"
#include "IObsWilsonBuilder.h"
#include "IObsCoreAPI.h"
#include "IWilsonFreezer.h"
#include "IObsQCDProxy.h"
#include "ObsWilsonHelper.h"

/**
 * @struct ObservablePortsConfig
 * @brief Dependency container (“ports”) for observable computations.
 *
 * This structure groups together all service-like interfaces required by
 * the observable layer to:
 *
 * - build Wilson coefficient blocks on demand,
 * - access input parameters (SM, flavor, etc.),
 * - query QCD running quantities,
 * - check global runtime flags (e.g. MARTY backend),
 * - freeze/unfreeze Wilson coefficient blocks when switching observable needs.
 *
 * It plays the role of a lightweight dependency-injection container for
 * observable calculators, avoiding tight coupling to concrete implementations.
 *
 * Typical usage:
 * @code
 * ObservablePortsConfig ports(
 *     obs_wilson_builder,
 *     obs_param_proxy_sm,
 *     obs_param_proxy_flavor,
 *     obs_qcd_proxy,
 *     obs_use_marty_api,
 *     obs_wilson_freezer
 * );
 * @endcode
 *
 * @see IObsWilsonBuilder
 * @see IObsParameterProxy
 * @see IObsQCDProxy
 * @see IWilsonFreezer
 */
struct ObservablePortsConfig {
    /**
     * @brief Construct the ports configuration.
     *
     * @param iobswb        Observable-level Wilson builder.
     * @param iobspp_sm     Proxy for SM / Wilson / observable parameters.
     * @param iobspp_flav   Proxy for flavor-specific parameters.
     * @param iobs_qcdp     Proxy for QCD running quantities.
     * @param iobs_use_marty API exposing whether the MARTY backend is active.
     * @param iobs_wfreezer Wilson freezer to freeze/unfreeze groups.
     */
    ObservablePortsConfig(std::shared_ptr<IObsWilsonBuilder> iobswb,
                          std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm,
                          std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav,
                          std::shared_ptr<IObsQCDProxy> iobs_qcdp,
                          std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty,
                          std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer,
                          std::shared_ptr<ObsWilsonHelper> iobs_whelper = nullptr) :
        iobswb(std::move(iobswb)),
        iobspp_sm(std::move(iobspp_sm)),
        iobspp_flav(std::move(iobspp_flav)),
        iobs_qcdp(std::move(iobs_qcdp)),
        iobs_use_marty(std::move(iobs_use_marty)),
        iobs_wfreezer(std::move(iobs_wfreezer)),
        iobs_whelper(std::move(iobs_whelper)) {
        if (!this->iobs_whelper) {
            this->iobs_whelper = std::make_shared<ObsWilsonHelper>();
        }
    }

    /// Builder used to (re)build Wilson coefficient groups on demand.
    std::shared_ptr<IObsWilsonBuilder> iobswb;

    /// Parameter proxy for SM / Wilson / observable parameters.
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm;

    /// Parameter proxy for flavor-specific parameters.
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_flav;

    /// Proxy giving access to QCD running quantities (alpha_s, running masses, etc.).
    std::shared_ptr<IObsQCDProxy> iobs_qcdp;

    /// Core API flag indicating whether the MARTY backend is in use.
    std::shared_ptr<IObsCoreAPI<bool>> iobs_use_marty;

    /// Freezer allowing Wilson groups to be frozen/unfrozen at the observable level.
    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer;

    /// Per-manager Wilson lifecycle helper. Must not be shared accidentally across independent managers.
    std::shared_ptr<ObsWilsonHelper> iobs_whelper;

};

#endif