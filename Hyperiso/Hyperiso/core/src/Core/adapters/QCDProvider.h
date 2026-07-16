#ifndef QCDADAPTER_H
#define QCDADAPTER_H

#include "IDataProvider.h"
#include "IQCDProvider.h"
#include "QCDHelper.h"
#include "Configs.h"

/**
 * @file QCDProvider.h
 * @brief High-level access to QCD quantities like α_s and running quark masses.
 *
 * This header declares the QCDProvider class, which implements both:
 *  - IDataProvider<QCDProvider> to fit into the generic data provider
 *    interface, and
 *  - IQCDProvider to expose access to QCDConstants.
 *
 * QCDProvider delegates its calculations to QCDHelper, which uses the
 * QCD-dependent block built from the SM inputs (see QCDHelper::Init()).
 */

/**
 * @class QCDProvider
 * @ingroup DataProvidersModule
 * @brief Provides strong coupling constant and MS-bar mass values at various energy scales.
 *
 * QCDProvider is a thin convenience wrapper around QCDHelper, offering a
 * callable interface for:
 *
 *  - computing the strong coupling \f$\alpha_s(\mu)\f$ from an AlphasConfig,
 *  - computing running MS-bar quark masses at scale \f$\mu\f$ from a MassConfig,
 *  - accessing the underlying QCDConstants used in these calculations.
 *
 * Typical usage:
 * @code
 * QCDProvider qcd;
 * AlphasConfig a_cfg{ .scale = 91.1876, .m_b_type = MassType::POLE, .m_t_type = MassType::POLE };
 * MassConfig   m_cfg{ .pdg_id = 5, .scale = 4.2, .m_b_type = MassType::POLE, .m_t_type = MassType::POLE };
 *
 * double alpha_s_mZ = qcd(a_cfg);
 * double mb_MSbar   = qcd(m_cfg);
 * auto*  consts     = qcd.get_constants();
 * @endcode
 */
class QCDProvider : public IDataProvider<QCDProvider>, public IQCDProvider {
public:
    /**
     * @brief Computes the strong coupling constant \f$\alpha_s\f$ at a given scale.
     *
     * Delegates to QCDHelper::alpha_s, using the configuration provided
     * (renormalization scale and mass-scheme choices for b and t quarks).
     *
     * @param config AlphasConfig object containing scale and mass type settings.
     * @return Value of \f$\alpha_s\f$ at the specified scale.
     */
    double operator()(AlphasConfig);

    /**
     * @brief Computes the MS-bar running mass at a given scale for a quark.
     *
     * Delegates to QCDHelper::msbar_mass using the supplied MassConfig
     * (PDG ID, renormalization scale, and mass schemes).
     *
     * @param config MassConfig object containing PDG ID and scale information.
     * @return The MS-bar mass at the specified scale.
     */
    double operator()(MassConfig);

    /**
     * @brief Retrieves the constants related to QCD calculations.
     *
     * Provides access to the QCDConstants structure used internally by
     * QCDHelper (Casimir operators, beta-function coefficients, etc.).
     *
     * @return Pointer to QCDConstants structure.
     */
    QCDConstants* get_constants() override;
};


#endif // QCDADAPTER_H
