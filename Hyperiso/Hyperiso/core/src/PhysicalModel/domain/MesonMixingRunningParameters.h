#ifndef __MESONMIXINGRUNNINGPARAMETERS_H__
#define __MESONMIXINGRUNNINGPARAMETERS_H__

#include <array>
#include "Math.h"

struct MesonMixingRunningParameters {
    static constexpr int n_coefs {8};
    static constexpr int n_pows {5};

    static constexpr std::array<double, n_pows> ai = {0.2609, 0.1304, -1.0435, -0.6315, 0.7184};
    static constexpr std::array<double, n_pows> bi = {0.2400, 0.1200, -0.6900, -0.5810, 0.6610};

    // For B0 mixing
    static constexpr double a0_V_5 = 1.;
    static constexpr double a1_V_5 = 1.6273;
    static constexpr double b_V_5 = -1.6273;

    static const std::array<std::array<std::array<double, 2>, 2>, 2> a0_LR_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> a1_LR_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> b_LR_5;

    static const std::array<std::array<std::array<double, 2>, 2>, 2> a0_S_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> a1_S_5;
    static const std::array<std::array<std::array<double, 2>, 2>, 2> b_S_5;

    // For K0/D0 mixing
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

	static inline std::array<complex_t, n_coefs> SUSY_to_BMU(const std::array<complex_t, n_coefs>& Ci_SUSY) {
		std::array<complex_t, n_coefs> Ci_BMU;
        for (size_t j = 0; j < n_coefs; j++) {
            for (size_t k = 0; k < n_coefs; k++) {
                Ci_BMU[j] += SUSY_to_BMU_superiso[j][k] * Ci_SUSY[k];
            }
        }
        return Ci_BMU;
	}

	static const std::array<std::array<double, n_coefs>, n_coefs> SUSY_to_BMU_superiso;
    static const std::array<std::array<double, n_coefs>, n_coefs> SUSY_to_BMU_fierz;
};

#endif // __MESONMIXINGRUNNINGPARAMETERS_H__
