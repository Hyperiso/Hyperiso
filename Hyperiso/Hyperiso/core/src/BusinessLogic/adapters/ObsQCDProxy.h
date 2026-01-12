#ifndef OBS_QCD_PROXY_H
#define OBS_QCD_PROXY_H

#include "IObsQCDProxy.h"
#include "QCDProvider.h"

/**
 * @file ObsQCDProxy.h
 * @brief Default observable-layer implementation of @ref IObsQCDProxy.
 *
 * This proxy delegates computations to @ref QCDProvider and exposes a stable,
 * simple API to the observable/business layer.
 *
 * @see IObsQCDProxy
 * @see QCDProvider
 */
class ObsQCDProxy : public IObsQCDProxy  {
public:
    /**
     * @brief Computes / retrieves alpha_s according to @ref AlphasConfig.
     */
    double operator()(AlphasConfig config) override;

    /**
     * @brief Computes / retrieves a (running) mass according to @ref MassConfig.
     */
    double operator()(MassConfig config) override;

    /**
     * @brief Returns a pointer to the internal QCD constants used by the provider.
     */
    QCDConstants *get_constants() override;

private:
    QCDProvider qcd_;
};

#endif // OBS_QCD_PROXY_H