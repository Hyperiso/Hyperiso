#ifndef MARTY_MASS_RUNNER_H
#define MARTY_MASS_RUNNER_H

#include "Parameters.h"

/**
 * @file MartyMassRunner.h
 * @brief Declares a lightweight wrapper around QCDHelper mass running.
 *
 * This header defines ::MartyMassRunner, a simple helper that exposes
 * running quark masses at a given matching scale to MARTY-generated
 * code in a convenient way.
 */

/**
 * @class MartyMassRunner
 * @ingroup DataProvidersModule
 * @brief Provides MS-bar running quark masses at a chosen scale.
 *
 * MartyMassRunner is a thin façade around ::QCDHelper::msbar_mass.
 * It stores a matching scale \f$Q_{\text{match}}\f$ and offers
 * direct getters for charm, bottom and top masses at that scale.
 */
class MartyMassRunner {
    /// Matching scale at which masses are evaluated.
    double Q_match;

public:
    /**
     * @brief Constructs the runner with a given matching scale.
     * @param Q_match Matching scale in GeV.
     */
    explicit MartyMassRunner(double Q_match) : Q_match(Q_match) {}

    /**
     * @brief Returns the \f$\overline{\text{MS}}\f$ top mass at @p Q_match.
     * @return Running top mass \f$m_t(Q_{\text{match}})\f$.
     */
    double get_mt() {
        return QCDHelper::msbar_mass(6, Q_match);
    }

    /**
     * @brief Returns the \f$\overline{\text{MS}}\f$ bottom mass at @p Q_match.
     * @return Running bottom mass \f$m_b(Q_{\text{match}})\f$.
     */
    double get_mb() {
        return QCDHelper::msbar_mass(5, Q_match);
    }

    /**
     * @brief Returns the \f$\overline{\text{MS}}\f$ charm mass at @p Q_match.
     * @return Running charm mass \f$m_c(Q_{\text{match}})\f$.
     */
    double get_mc() {
        return QCDHelper::msbar_mass(4, Q_match);
    }
};

#endif // MARTY_MASS_RUNNER_H
