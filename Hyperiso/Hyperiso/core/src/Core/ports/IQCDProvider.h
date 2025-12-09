#ifndef IQCDPROVIDER_H
#define IQCDPROVIDER_H

#include "AbstractConfig.h"
#include "QCDHelper.h"

/**
 * @file IQCDProvider.h
 * @brief Abstract interface for accessing QCD constants.
 *
 * This header declares the IQCDProvider interface, which exposes an abstract
 * way to retrieve a QCDConstants structure. It is meant to be implemented
 * by concrete data providers such as QCDProvider, which perform actual QCD
 * calculations and want to make their internal constants available.
 *
 * ## Related Classes
 * - @ref QCDConstants
 * - @ref QCDHelper
 * - @ref QCDProvider
 */

/**
 * @class IQCDProvider
 * @ingroup DataProvidersModule
 * @brief Abstract interface for accessing QCD constants.
 *
 * Implementations of this interface provide read access to the internal
 * QCDConstants used to compute:
 *  - the strong coupling \f$\alpha_s\f$,
 *  - running quark masses,
 *  - threshold matching coefficients, etc.
 *
 * It is intentionally small and focused, and is usually combined with a
 * more feature-complete data provider (e.g. QCDProvider).
 */
class IQCDProvider {
public:
    /// Virtual destructor for interface.
    virtual ~IQCDProvider() = default;

    /**
     * @brief Retrieves the QCD constants used by the implementing provider.
     *
     * The returned pointer is expected to remain valid for the lifetime of
     * the implementing object (typically it points to QCDHelper::constants).
     *
     * @return Pointer to a QCDConstants structure.
     */
    virtual QCDConstants* get_constants() = 0;
};

#endif // IQCDPROVIDER_H
