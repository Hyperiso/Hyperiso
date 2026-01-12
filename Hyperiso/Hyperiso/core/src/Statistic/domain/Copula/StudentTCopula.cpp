#include "StudentTCopula.h"

StudentTCopula::StudentTCopula(unsigned int seed, RealMatrix R, int nu) : GenericCopula(seed) {
    if (nu < 2)
        throw std::invalid_argument("Number of DoF should be at least 2 for finite variance.");

    this->nu = nu;
    this->R = nearest_psd(R);
    this->L = cholesky_L(R);
    this->logdet = R.slogdet().logdet;
    this->R_inv = R.inv();
}

std::vector<Vector> StudentTCopula::sample_u(std::size_t n) {
    std::size_t d = R.cols();
    std::vector<std::vector<double>> U;

    RealMatrix z (d, 1);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < d; j++) {
            z.at(j, 0) = gsl_ran_ugaussian(eng_);   // z follows MN(0, 1)
        }

        z = L * z; // z follows MN(0, R)
        double w = gsl_ran_chisq(eng_, nu);
        z /= std::sqrt(w / nu); // z follows Mt(R, nu)

        Vector u (d, 0.0);
        for (std::size_t j = 0; j < d; j++) {
            u[j] = std::clamp(gsl_cdf_tdist_P(z.at(j, 0), nu), CLIP_U, 1 - CLIP_U);   // u follows C(R, nu)
        }
        U.emplace_back(std::move(u));
    }

    return U;
}

double StudentTCopula::log_density(Vector u) {
    std::size_t d = u.size();
    RealMatrix z (d, 1);

    double log_t1 {0.0};
    for (size_t i = 0; i < d; i++) {
        z.at(i, 0) = gsl_cdf_tdist_Pinv(u[i], nu);
        log_t1 += std::log(gsl_ran_tdist_pdf(z.at(i, 0), nu));
    }

    double quad = (z.transpose() * R_inv * z).at(0, 0);
    double log_td = gsl_sf_lngamma((nu + d) / 2) - gsl_sf_lngamma(nu / 2)
                    - (d / 2) * std::log(nu * PI) - 0.5 * logdet
                    - (nu + d) / 2 * gsl_sf_log_1plusx(quad / nu);

    return log_t1 + log_td;
}
