#ifndef IOBSWILSONBUILDER_H
#define IOBSWILSONBUILDER_H

#include <memory>

#include "AbstractConfig.h"

class IObsWilsonProxy;

/**
 * @file IObsWilsonBuilder.h
 * @brief Interface for building Wilson coefficient services in the observable layer.
 *
 * This is the observable-side counterpart to the core Wilson builder system.
 * Implementations:
 *  - accept a generic @ref AbstractConfig (runtime config object),
 *  - build (or extend) the Wilson coefficient infrastructure,
 *  - provide an @ref IObsWilsonProxy for coefficient queries.
 *
 * @see ObsWilsonBuilder
 * @see IObsWilsonProxy
 */

class IObsWilsonBuilder {
public: 
    virtual ~IObsWilsonBuilder() = default;	

    /**
     * @brief Build or extend the Wilson coefficient infrastructure using an observable config.
     * @param config Configuration object (typically a @ref WilsonBuildConfig).
     */
    virtual void build(std::shared_ptr<AbstractConfig>) = 0;

    /**
     * @brief Return a proxy used by observables to access Wilson coefficients.
     * @return An observable Wilson proxy.
     */
    virtual std::shared_ptr<IObsWilsonProxy> get_proxy() = 0;
};

#endif // IOBSWILSONBUILDER_H
