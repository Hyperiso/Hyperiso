#if !defined(HYPERISO_FORMFACTORS_H)
#define HYPERISO_FORMFACTORS_H

/*
*   Numerical constants for flavor observables 
*/ 

constexpr double GAMMA = 0.5772156649015328;

// Digamma function for integer inputs
static double psi(int n) {
    double sum = 0;
    for (int k = 1; k < n; ++k) {
        sum += 1 / k;
    }
    return sum - GAMMA;
}

namespace FFInput {

    // B -> K* gamma isospin asymmetry

    // Gegenbauer expansion coefficients at mu = 1 GeV for the Light-cone wavefunctions [ref missing]
    constexpr double a_1_perp = 0.04;
    constexpr double a_2_perp = 0.10;
    constexpr double a_1_par = 0.06;
    constexpr double a_2_par = 0.16; 

    // Anomalous dimensions for the Gegenbauer expansion coefficients [hep-ph/9802299]
    static double gamma_n_perp(int n, int Nc) {
        double C_F = (Nc * Nc - 1) / (2 * Nc);
        return 4 * C_F * (psi(n + 1) + GAMMA - .75 + 1 / (n + 1));
    }

    static double gamma_n_par(int n, int Nc) {
        double C_F = (Nc * Nc - 1) / (2 * Nc);
        return 4 * C_F * (psi(n + 2) + GAMMA - .75 - 1 /(2 * (n + 1) * (n + 2)));
    }

    // Other Gegenbauer momenta parameters [ref missing]
    constexpr double zeta_3_V = 0.032;
    constexpr double zeta_3_A = 0.013;
    constexpr double omega_10_A = -2.1;
    constexpr double delta_tilde_p = 0.16;
    constexpr double delta_tilde_m = -0.16;

    // Parameter for the first inverse moment of the B meson distribution [ref missing]
    constexpr double lambda_B = 0.46; 

    // T1 formfactor (0.282 LCSR only, 0.312 LCSR + Lattice, 0.35 Beneke et al. [hep-ph/0106067])
    constexpr double T1_B_Kstar = 0.312;

}

#endif // HYPERISO_FORMFACTORS_H