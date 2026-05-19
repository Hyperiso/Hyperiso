#ifndef CHISQUAREDLIKELIHOOD_H
#define CHISQUAREDLIKELIHOOD_H

#include "BaseLikelihood.h"

#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * Likelihood chi2 global sans nuisances explicites.
 *
 * theta == p, eta est vide.
 *
 * nll(p) = 1/2 * (f(p, empty_eta) - O_exp)^T C^{-1} (f(p, empty_eta) - O_exp)
 *
 * La matrice C^{-1} est typiquement l'inverse de la covariance MC des prédictions
 * autour du point SM / central.
 */
class ChiSquaredLikelihood final : public BaseLikelihood {
public:
    ChiSquaredLikelihood(
        const ModelFn& model,
        std::shared_ptr<LikelihoodContext> ctx,
        std::size_t p_dim,
        RealMatrix covariance_inv
    );

    double nll(const std::vector<double>& theta) const override;
    std::size_t dim() const override;

    RealMatrix observable_curvature(
        const std::vector<double>& r
    ) const override;

    RealMatrix nuisance_curvature(
        const std::vector<double>& eta
    ) const override;

private:
    RealMatrix covariance_inv_;
};

#endif // CHISQUAREDLIKELIHOOD_H
