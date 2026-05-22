#include "ChiSquaredLikelihood.h"

#include <iostream>

ChiSquaredLikelihood::ChiSquaredLikelihood(
    const ModelFn& model,
    std::shared_ptr<LikelihoodContext> ctx,
    std::size_t p_dim,
    RealMatrix covariance_inv
)
    : BaseLikelihood(model, std::move(ctx), p_dim),
      covariance_inv_(std::move(covariance_inv))
{}

double ChiSquaredLikelihood::nll(const std::vector<double>& theta) const {
    if (theta.size() != p_dim) {
        throw std::invalid_argument("ChiSquaredLikelihood::nll: theta size must equal p_dim");
    }

    try {
        const std::vector<double> eta_empty;
        std::vector<double> pred = model(theta, eta_empty);

        if (pred.size() != ctx->exp_obs_values.size()) {
            throw std::runtime_error("ChiSquaredLikelihood::nll: prediction/observation size mismatch");
        }
        if (covariance_inv_.rows() != pred.size() || covariance_inv_.cols() != pred.size()) {
            throw std::runtime_error("ChiSquaredLikelihood::nll: covariance inverse dimension mismatch");
        }

        std::vector<double> r(pred.size(), 0.0);
        for (std::size_t i = 0; i < pred.size(); ++i) {
            r[i] = pred[i] - ctx->exp_obs_values[i];
            if (!std::isfinite(r[i])) {
                return 1e100;
            }
        }

        double q = 0.0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            double row = 0.0;
            for (std::size_t j = 0; j < r.size(); ++j) {
                const double cij = covariance_inv_.at(i, j);
                if (!std::isfinite(cij)) {
                    return 1e100;
                }
                row += cij * r[j];
            }
            q += r[i] * row;
        }

        const double out = 0.5 * q; //Niels avait tort
        // const double out = q; //Niels à dit pas de 1/2, c'est sa responsabilité.
        return std::isfinite(out) ? out : 1e100;
    } catch (const std::exception& e) {
        std::cout << "[CHI2DBG] nll exception: " << e.what() << "\n";
        return 1e100;
    } catch (...) {
        std::cout << "[CHI2DBG] nll unknown exception\n";
        return 1e100;
    }
}

std::size_t ChiSquaredLikelihood::dim() const {
    return p_dim;
}

RealMatrix ChiSquaredLikelihood::observable_curvature(
    const std::vector<double>&
) const {
    return covariance_inv_;
}

RealMatrix ChiSquaredLikelihood::nuisance_curvature(
    const std::vector<double>& eta
) const {
    if (!eta.empty()) {
        throw std::invalid_argument("ChiSquaredLikelihood::nuisance_curvature: eta must be empty");
    }
    return RealMatrix(0, 0);
}
