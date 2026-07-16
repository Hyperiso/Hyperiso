#ifndef MESONMIXINGRUNNINGPARAMETERS_H
#define MESONMIXINGRUNNINGPARAMETERS_H

#include <array>
#include <cstddef>

#include "Math.h"

/**
 * @file MesonMixingRunningParameters.h
 * @brief Numerical constants and basis-change utilities for meson-mixing running matrices.
 *
 * This struct groups all constants used to build the running evolution matrices
 * relevant to meson mixing (B0, K0, D0), including:
 *  - evolution exponents (ai, bi),
 *  - coefficient matrices for LL/LR/S sectors at nf=5 and nf=4 steps,
 *  - basis change matrices (SUSY <-> BMU),
 *  - a helper to apply a basis change to a coefficient vector.
 *
 * The actual values are defined in the corresponding .cpp to keep this header lightweight.
 */
struct MesonMixingRunningParameters {
    /// Number of Wilson coefficients for the mixing operator basis handled here.
    static constexpr int n_coefs {8};

    /// Number of exponent entries stored in ai / bi (for eta power vectors).
    static constexpr int n_pows {5};

    /// Exponents for nf=5 evolution (eta_5^ai).
    static constexpr std::array<double, n_pows> ai = {0.2609, 0.1304, -1.0435, -0.6315, 0.7184};

    /// Exponents for nf=4 evolution (eta_4^bi).
    static constexpr std::array<double, n_pows> bi = {0.2400, 0.1200, -0.9600, -0.5810, 0.6610};

    // -------------------------------------------------------------------------
    // B0 mixing (nf=5 evolution) : V, LR, S sectors
    // -------------------------------------------------------------------------
    
    static constexpr double a0_V_5 = 1.;
    static constexpr double a1_V_5 = 1.6273;
    static constexpr double b_V_5 = -1.6273;

    static const std::array<std::array<std::array<double, 2>, 2>, 2> a0_LR_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> a1_LR_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> b_LR_5;

    static const std::array<std::array<std::array<double, 2>, 2>, 2> a0_S_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> a1_S_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> b_S_5;

    // -------------------------------------------------------------------------
    // K0 / D0 mixing (nf=4 evolution) : V, LR, S sectors
    // -------------------------------------------------------------------------

    static constexpr double a0_V_4 = 1.;
    static constexpr double a1_V_4 = 1.7917;
    static constexpr double b_V_4 = -0.1644;
    static constexpr double c_V_4 = -1.6273;

    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> a0_LR_4;
    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> a1_LR_4;
    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> b_LR_4;
    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> c_LR_4;

    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> a0_S_4;
    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> a1_S_4;
    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> b_S_4;
    static const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> c_S_4;

    // -------------------------------------------------------------------------
    // Basis changes
    // -------------------------------------------------------------------------

    /**
     * @brief Applies a basis change Ci_out = P * Ci_in.
     *
     * @param Ci_in Vector of coefficients in the input basis.
     * @param P     Change-of-basis matrix.
     * @return      Vector of coefficients in the output basis.
     */
	static inline std::array<complex_t, n_coefs> change_basis(const std::array<complex_t, n_coefs>& Ci_in, const std::array<std::array<double, n_coefs>, n_coefs>& P) {
		std::array<complex_t, n_coefs> Ci_out;
        for (size_t j = 0; j < n_coefs; j++) {
            for (size_t k = 0; k < n_coefs; k++) {
                Ci_out[j] += P[j][k] * Ci_in[k];
            }
        }
        return Ci_out;
	}

    /// Change of basis matrix from SUSY basis to BMU basis.
	static const std::array<std::array<double, n_coefs>, n_coefs> SUSY_to_BMU;

    /// Change of basis matrix from BMU basis to SUSY basis.
    static const std::array<std::array<double, n_coefs>, n_coefs> BMU_to_SUSY;
};

#endif // MESONMIXINGRUNNINGPARAMETERS_H
