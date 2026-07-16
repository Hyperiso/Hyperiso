#ifndef QCDHELPER_H
#define QCDHELPER_H

#include <array>
#include <string>

#include "Parameters.h"
#include "Math.h"
#include "DependentBlockManager.h"

/**
 * @file QCDHelper.h
 * @brief Core QCD utilities for α_s running and quark mass evolution.
 *
 * This header declares:
 *  - QCDConstants: a collection of static QCD coefficients (Casimirs, beta
 *    and gamma function coefficients).
 *  - QCDHelper: the numerical engine computing:
 *      * the running strong coupling α_s(μ),
 *      * running MS-bar quark masses,
 *      * various QCD scales (Λ_QCD for different nf),
 *      * pole and threshold-related masses.
 *
 * QCDHelper also defines and populates the derived "QCD" block via a
 * DependentBlock using DependentBlockManager and SM parameters.
 */

 /**
 * @struct QCDConstants
 * @brief Collection of QCD-specific constants and coefficients.
 *
 * Contains:
 *  - Color factors for SU(3): N_c, C_F, C_A,
 *  - β-function coefficients (up to three loops) for nf = 1..6,
 *  - γ-function (anomalous dimension) coefficients for the quark mass
 *    running (up to three loops).
 *
 * These constants are exposed through QCDHelper::constants and used in
 * the explicit α_s and mass-running formulae.
 */
struct QCDConstants {
    /// Number of colors.
    static constexpr int Nc = 3;
    /// Fundamental Casimir C_F = (N_c^2 - 1) / (2 N_c).
    static constexpr double C_F = (Nc * Nc - 1.) / (2. * Nc);
    /// Adjoint Casimir C_A = N_c.
    static constexpr double C_A = Nc;

    /**
     * @brief β-function coefficients β_0, β_1, β_2 for nf = 1..6.
     *
     * Indexed as beta[nf-1][i] with i = 0,1,2.
     */
    static constexpr std::array<std::array<double, 3>, 6> beta {{{31. / 3, 134. / 3, 2309.8}, 
                                                                 {29. / 3, 115. / 3, 1786.7}, 
                                                                 {9      , 32      , 1287.6}, 
                                                                 {25. / 3, 77. / 3 , 812.7 }, 
                                                                 {23. / 3, 58. / 3 , 361.81}, 
                                                                 {7      , 13      , -65   }}};
          
    /**
     * @brief γ-function coefficients γ_0, γ_1, γ_2 for nf = 1..6.
     *
     * Indexed as gamma[nf-1][i] with i = 0,1,2.
     */
    static constexpr std::array<std::array<double, 3>, 6> gamma {{{2, 8.1388, 34.408},
                                                                  {2, 7.8611, 29.678}, 
                                                                  {2, 7.5833, 24.840}, 
                                                                  {2, 7.3055, 19.894}, 
                                                                  {2, 7.0277, 14.839}, 
                                                                  {2, 6.75  , 9.6773}}};
};

/**
 * @class QCDHelper
 * @brief Static utility class for QCD running and derived QCD block construction.
 *
 * QCDHelper:
 *  - builds the "QCD" dependent block from SM inputs via Init(),
 *  - provides α_s(μ) at arbitrary scales and for different mass schemes,
 *  - provides running MS-bar masses for quarks,
 *  - abstracts the matching across flavor thresholds (nf changes),
 *  - exposes the internal QCDConstants via a static pointer.
 *
 * It is intended as the low-level engine behind higher-level providers like
 * QCDProvider.
 */
class QCDHelper {
private:
    /// Internal mass scheme used for b-quark in threshold/matching (default pole).
    static inline MassType m_b_type {MassType::POLE};
    /// Internal mass scheme used for t-quark in threshold/matching (default pole).
    static inline MassType m_t_type {MassType::POLE};

    /// Returns Λ_QCD( nf ) chosen such that α_s(μ) matches the input configuration.
    static double get_lambda(double mu, MassType mass_b_type, MassType mass_t_type);

    /// Returns ordered quark-mass thresholds to decide active flavors nf.
    static std::vector<double> getOrderedMasses(MassType mass_b_type, MassType mass_t_type);

    /// Finds Λ that reproduces a target α_s at scale Q for a given nf.
    static double match_lambda(double target_alpha, double Q, int nf);

    /// Explicit multi-loop expression for α_s(μ; Λ, nf).
    static double alpha_s_explicit(double mu, double lambda, int nf);

    /// Runs a quark mass from Q_i to Q_f with nf active flavors.
    static double runMass(double mass, double Q_i, double Q_f, int nf, MassType m_b_type, MassType m_t_type);

    /// Helper factor for running mass solution (R(α, nf)).
    static double R(double alpha, int nf);

    /// Computes charm pole mass from Λ_4 and SM inputs.
    static double calc_mc_pole(double lambda_4);
    /// Computes charm pole mass from Λ_4 and SM inputs at one loop.
    static double calc_mc_pole_one_loop(double lambda_4);
    /// Computes bottom pole mass from Λ_5 and SM inputs.
    static double calc_mb_pole(double lambda_5);
    /// Computes bottom pole mass from Λ_5 and SM inputs at one loop.
    static double calc_mb_pole_one_loop(double lambda_5);
    /// Computes kinematic bottom mass from m_b(m_b) and Λ_3.
    static double calc_mb_kinematic(double mb_mb, double lambda_3);
    /// Computes 1S bottom mass from Λ_4 and m_b^pole.
    static double calc_mb_1S(double lambda_4, double mb_pole);
    /// Computes m_t(m_t) from Λ_6 and Λ_5 (two-step procedure).
    static double calc_mt_mt(double lambda6_mt_pole, double lambda_5);

public:
    /// Static storage for QCD constants.
    static inline QCDConstants _static_constants{};
    /// Pointer to QCD constants (to conform with IQCDProvider interface).
    static inline QCDConstants* constants = &_static_constants;

    /**
     * @brief Initializes the QCD dependent block.
     *
     * Creates and registers a DependentBlock named "QCD" using
     * DependentBlockManager, taking SM inputs (SMINPUTS, MASS) as sources.
     * This block stores quantities such as:
     *  - Λ_3, Λ_4, Λ_5, Λ_6,
     *  - m_b in various schemes,
     *  - m_t(m_t), etc.
     *
     * Must be called once during the SM initialization sequence
     * (see SMModelStrategy::postInitialization).
     */
    static void Init();

    /**
     * @brief Computes the strong coupling constant α_s at scale μ.
     *
     * Uses the internally constructed QCD block (Λ_n) and the SM quark masses
     * to determine the appropriate nf and perform matching across thresholds.
     *
     * @param mu          Renormalization scale μ.
     * @param mass_b_type Mass scheme for bottom (pole or MS-bar) in thresholds.
     * @param mass_t_type Mass scheme for top (pole or MS-bar) in thresholds.
     * @return α_s(μ).
     */
    static double alpha_s(double mu, MassType mass_b_type = MassType::POLE, MassType mass_t_type = MassType::POLE);

    /**
     * @brief Computes the MS-bar running mass of a quark at scale μ.
     *
     * The starting point (initial mass and scale) is inferred from SM / QCD
     * blocks, then the mass is evolved to μ using the appropriate nf and the
     * R(α, nf) factor.
     *
     * @param pdg_code    PDG code of the quark (1..6).
     * @param mu          Target renormalization scale.
     * @param mass_b_type Mass scheme for bottom thresholds.
     * @param mass_t_type Mass scheme for top thresholds.
     * @return m_q(μ) in MS-bar scheme.
     */
    static double msbar_mass(int pdg_code, double mu, MassType mass_b_type = MassType::POLE, MassType mass_t_type = MassType::POLE);

    /**
     * @brief Returns the number of active flavors n_f at scale μ.
     *
     * Uses ordered thresholds derived from the chosen mass schemes
     * (m_b_type, m_t_type).
     *
     * @param mu          Renormalization scale.
     * @param mass_b_type Bottom mass scheme.
     * @param mass_t_type Top   mass scheme.
     * @return Number of active flavors n_f.
     */
    static int get_nf(double mu, MassType mass_b_type = MassType::POLE, MassType mass_t_type = MassType::POLE);
};

#endif // QCDHELPER_H
