#include "GaussianCopula.h"

GaussianCopula::GaussianCopula(unsigned int seed, RealMatrix R) : GenericCopula(seed) {
    R = nearest_psd(R);
    L = cholesky_L(R);
    logdet = R.slogdet().logdet;
    R_inv = R.inv();
}

std::vector<Vector> GaussianCopula::sample_u(std::size_t n) {
    std::vector<std::vector<double>> U;

    for (std::size_t i = 0; i < n; i++) {
        U.emplace_back(sample_u());
    }

    return U;
}

Vector GaussianCopula::sample_u() {
    std::size_t d = L.cols();

    RealMatrix z (d, 1);
    for (std::size_t j = 0; j < d; j++) {
        z.at(j, 0) = gsl_ran_ugaussian(eng_);   // z follows MN(0, 1)
    }

    z = L * z; // z follows MN(0, R)
    Vector u (d, 0.0);
    for (std::size_t j = 0; j < d; j++) {
        u[j] = std::clamp(gsl_cdf_ugaussian_P(z.at(j, 0)), CLIP_U, 1 - CLIP_U);   // u follows C(R)
    }

    return u;
}

double GaussianCopula::log_density(Vector u) {
    std::size_t d = u.size();
    RealMatrix z (d, 1);

    for (size_t i = 0; i < d; i++) {
        z.at(i, 0) = gsl_cdf_ugaussian_Pinv(std::clamp(u[i], 1e-13, 1. - 1e-13));
    }

    // return -0.5 * logdet - 0.5 * (z.transpose() * (R_inv - eye(d)) * z).at(0, 0);
    return -0.5 * (z.transpose() * (R_inv - eye(d)) * z).at(0, 0);
}
