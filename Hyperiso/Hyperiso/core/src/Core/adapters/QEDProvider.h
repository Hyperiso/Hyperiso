#ifndef QEDPROVIDER_H
#define QEDPROVIDER_H

#include "IDataProvider.h"
#include "EWHelper.h"
#include "Configs.h"

/**
 * @file QEDProvider.h
 * @brief High-level access to QED/electroweak quantities like the electromagnetic coupling.
 *
 * This header declares the QEDProvider class, which implements
 * IDataProvider<QEDProvider> to fit into the generic data provider interface.
 *
 * QEDProvider delegates its calculations to EWHelper, which builds and uses
 * the EW-dependent block from SM, QCD, and Wilson-scale inputs
 * (see EWHelper::Init()).
 */

/**
 * @class QEDProvider
 * @ingroup DataProvidersModule
 * @brief Provides electromagnetic coupling values at requested energy scales.
 *
 * QEDProvider is a thin convenience wrapper around EWHelper, offering a
 * callable interface for computing the electromagnetic coupling
 * \f$\alpha_{\text{em}}(\mu)\f$ from an AlphasConfig.
 *
 * @note EWHelper::alpha_em does not implement generic running and throws
 * std::logic_error instead of returning a placeholder value. Code requiring
 * the precomputed special values should read them from
 * the EW block populated by EWHelper::Init().
 *
 * Typical usage:
 * @code
 * QEDProvider qed;
 * AlphasConfig a_cfg{ .scale = 91.1876 };
 *
 * double alpha_em = qed(a_cfg);
 * @endcode
 */
class QEDProvider : public IDataProvider<QEDProvider> {
public:
    /**
     * @brief Computes the electromagnetic coupling \f$\alpha_{\text{em}}\f$ at a given scale.
     *
     * Delegates to EWHelper::alpha_em using the renormalization scale stored
     * in the supplied configuration.
     *
     * @param config AlphasConfig object containing the target scale.
     * @return Value of \f$\alpha_{\text{em}}\f$ at the specified scale, as returned by EWHelper.
     */
    double operator()(AlphasConfig);
};

#endif // QEDPROVIDER_H
