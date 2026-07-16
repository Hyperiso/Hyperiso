#ifndef IOBS_QCD_PROXY_H
#define IOBS_QCD_PROXY_H

#include "QCDProvider.h"

/**
 * @file IObsQCDProxy.h
 * @brief Interface for accessing QCD-derived quantities from the observable layer.
 *
 * The observable/business layer sometimes needs QCD "services" without directly
 * coupling to low-level providers. This interface defines the minimal contract
 * to retrieve:
 *  - alpha_s(μ) or related quantities (via @ref AlphasConfig),
 *  - quark masses or running masses (via @ref MassConfig),
 *  - the underlying QCD constants container (if needed by advanced routines).
 *
 * Typical usage:
 * @code
 *   std::shared_ptr<IObsQCDProxy> qcd = ...;
 *   double as_mw = (*qcd)(AlphasConfig{...});
 *   double mb    = (*qcd)(MassConfig{...});
 * @endcode
 *
 * Concrete implementations usually wrap @ref QCDProvider.
 *
 * @see ObsQCDProxy
 * @see QCDProvider
 * @see AlphasConfig
 * @see MassConfig
 */
class IObsQCDProxy {
public:
    virtual ~IObsQCDProxy() = default;

    /**
     * @brief Computes / retrieves alpha_s according to the given configuration.
     * @param config Configuration describing scale, scheme, thresholds, etc.
     * @return Requested strong coupling value.
     */
    virtual double operator()(AlphasConfig config) = 0;

    /**
     * @brief Computes / retrieves a (running) quark mass according to the given configuration.
     * @param config Configuration describing flavor, scale, scheme, etc.
     * @return Requested mass value.
     */
    virtual double operator()(MassConfig config) = 0;

    /**
     * @brief Returns a pointer to the internal QCD constants.
     *
     * This is intended for advanced use cases. Most clients should prefer calling
     * the operator() overloads.
     *
     * @return Pointer to an internal constants container (lifetime owned by implementation).
     */
    virtual QCDConstants* get_constants() = 0;
};

#endif // IOBS_QCD_PROXY_H
