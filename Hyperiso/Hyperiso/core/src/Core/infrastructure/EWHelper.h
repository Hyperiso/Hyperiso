#ifndef EWHELPER_H
#define EWHELPER_H

#include <array>
#include <string>

#include "Parameters.h"
#include "Math.h"
#include "analysis.h"
#include "QCDHelper.h"

/**
 * @file EWHelper.h
 * @brief Electroweak helper utilities and dependent-block initialization.
 *
 * This header declares:
 *  - ::QCDHelper–compatible routines for electroweak (EW) quantities,
 *  - the ::EWHelper class that creates a dependent "EW" block via
 *    DependentBlockManager and provides access to effective \f$\alpha_{\text{em}}\f$
 *    values at specific scales.
 *
 * All core formulas are based on hep-ph/0011135.
 */

// All calculations are taken from hep-ph/0011135

/**
 * @class EWHelper
 * @ingroup SpectrumCalculationModule
 * @brief Helper class for electroweak quantities and EW-dependent blocks.
 *
 * EWHelper is responsible for:
 *  - declaring and initializing a dependent "EW" block through
 *    DependentBlockManager::addDependentBlock,
 *  - providing (currently stubbed) running of the electromagnetic coupling
 *    \f$\alpha_{\text{em}}(\mu)\f$,
 *  - implementing internal building blocks (delta-leptonic/hadronic,
 *    matching across thresholds, etc.) for alpha evolution.
 *
 * The "EW" block stores several special values of \f$\alpha_{\text{em}}\f$:
 *  - (1,1): \f$\alpha_{\text{em}}(m_Z)\f$,
 *  - (1,2): \f$\alpha_{\text{em}}(\mu_b)\f$,
 *  - (1,3): \f$\alpha_{\text{em}}(\mu_c)\f$,
 *  - (1,4): \f$\alpha_{\text{em}}(0) \approx 1/137.036\f$ (Thomson limit).
 *
 * These values are computed from the content of SMINPUTS, QCD, and Wilson
 * scale blocks (B_SCALE, D_SCALE).
 */
class EWHelper {
private:
    /**
     * @brief QCD-weighted integral used in EW hadronic contributions.
     *
     * @param k        Power index for the integrand.
     * @param mu_low   Lower integration scale.
     * @param mu_high  Upper integration scale.
     * @return Value of the integral \f$\int_{\mu_{\text{low}}^2}^{\mu_{\text{high}}^2} \alpha_s(\sqrt{s}) / s \, ds\f$.
     */
    static double I_s(int k, double mu_low, double mu_high);

    /**
     * @brief Leptonic contribution to the running of \f$\alpha_{\text{em}}\f$.
     *
     * @param mu_low   Lower scale.
     * @param mu_high  Upper scale.
     * @param alpha    Value of \f$\alpha_{\text{em}}\f$ used in the correction.
     * @return Leptonic contribution \f$\Delta_{\text{lept}}(\mu_{\text{low}}, \mu_{\text{high}})\f$.
     */
    static double delta_lept(double mu_low, double mu_high, double alpha);

    /**
     * @brief Perturbative (partonic) contribution to the running of \f$\alpha_{\text{em}}\f$.
     *
     * @param mu_low   Lower scale.
     * @param mu_high  Upper scale.
     * @param n_f      Number of active quark flavors.
     * @param alpha    Value of \f$\alpha_{\text{em}}\f$.
     * @return Partonic contribution \f$\Delta_{\text{part}}(\mu_{\text{low}}, \mu_{\text{high}})\f$.
     */
    static double delta_part(double mu_low, double mu_high, int n_f, double alpha);

    /**
     * @brief Hadronic contribution to the running of \f$\alpha_{\text{em}}\f$.
     *
     * @param mu_low   Lower scale.
     * @param mu_high  Upper scale.
     * @param n_f      Number of active quark flavors.
     * @param alpha    Value of \f$\alpha_{\text{em}}\f$.
     * @return Hadronic contribution \f$\Delta_{\text{had}}(\mu_{\text{low}}, \mu_{\text{high}})\f$.
     */
    static double delta_had(double mu_low, double mu_high, int n_f, double alpha);

    /**
     * @brief Heavy-quark vacuum polarization contribution.
     *
     * @param mu        Renormalization scale.
     * @param e_q       Electric charge of the heavy quark.
     * @param m_q_pole  Pole mass of the heavy quark.
     * @param n_l       Number of light flavors.
     * @param alpha     \f$\alpha_{\text{em}}\f$ used in the correction.
     * @return Heavy-quark contribution \f$\Pi_{\text{heavy}}\f$.
     */
    static double pi_heavy(double mu, double e_q, double m_q_pole, double n_l, double alpha);

    /**
     * @brief Mixed light–heavy contribution to the vacuum polarization.
     *
     * @param mu        Renormalization scale.
     * @param m_q_pole  Pole mass of the heavy quark.
     * @param n_l       Number of light flavors.
     * @return Mixed light–heavy contribution.
     */
    static double pi_light_heavy(double mu, double m_q_pole, double n_l);

    /**
     * @brief Computes \f$\alpha_{\text{em}}\f$ at \f$\mu_b\f$ with 5 active flavors.
     *
     * @param inv_alpha_m_Z  Inverse \f$\alpha_{\text{em}}\f$ at the Z pole.
     * @param m_Z            Z-boson mass.
     * @param mu_b           Bottom threshold scale.
     */
    static double alpha_5_mu_b(double inv_alpha_m_Z, double m_Z, double mu_b);

    /**
     * @brief Computes \f$\alpha_{\text{em}}\f$ at \f$\mu_b\f$ with 4 active flavors.
     *
     * @param inv_alpha_m_Z  Inverse \f$\alpha_{\text{em}}\f$ at the Z pole.
     * @param m_Z            Z-boson mass.
     * @param m_b_pole       Bottom pole mass.
     * @param mu_b           Bottom threshold scale.
     */
    static double alpha_4_mu_b(double inv_alpha_m_Z, double m_Z, double m_b_pole, double mu_b);

    /**
     * @brief Computes \f$\alpha_{\text{em}}\f$ at \f$\mu_c\f$ with 4 active flavors.
     *
     * @param inv_alpha_m_Z  Inverse \f$\alpha_{\text{em}}\f$ at the Z pole.
     * @param m_Z            Z-boson mass.
     * @param m_b_pole       Bottom pole mass.
     * @param mu_b           Bottom threshold scale.
     * @param mu_c           Charm threshold scale.
     */
    static double alpha_4_mu_c(double inv_alpha_m_Z, double m_Z, double m_b_pole, double mu_b, double mu_c);

public:
    /**
     * @brief Initializes the EW dependent block.
     *
     * This constructs the "EW" DependentBlock using DependentBlockManager::addDependentBlock,
     * with the following sources:
     *  - SM::SMINPUTS, SM::QCD,
     *  - WILSON::B_SCALE, WILSON::D_SCALE.
     *
     * The resulting block contains several fixed \f$\alpha_{\text{em}}\f$ values at:
     *  - \f$m_Z\f$,
     *  - \f$\mu_b\f$,
     *  - \f$\mu_c\f$,
     *  - and the Thomson limit (\f$\mu \to 0\f$).
     */
    static void Init();

    /**
     * @brief Returns \f$\alpha_{\text{em}}(\mu)\f$ at a given scale.
     *
     * Currently, the full running is not implemented, and this function
     * simply logs an error and returns 0.0. Code that needs \f$\alpha_{\text{em}}\f$
     * should instead use the precomputed special values in the "EW" block.
     *
     * @param mu Renormalization scale.
     * @return Value of \f$\alpha_{\text{em}}(\mu)\f$ (currently always 0.0 with an error log).
     */
    static double alpha_em(double mu);
};

#endif // EWHELPER_H
